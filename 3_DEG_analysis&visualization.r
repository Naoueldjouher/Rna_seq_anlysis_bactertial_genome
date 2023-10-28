library("GenomicFeatures")
library("lapply")
library("Rsamtools")
library("rtracklayer")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("ggplot2")
setwd("./Genome/Bacterial_rna_seq/Data/Output")




db2 <- makeTxDbFromGFF(file="./Genome/Bacterial_rna_seq/GCA_002173515.1_ASM217351v1_genomic.gff",
                      format="gff", organism="Komagataeibacter europaeus",
                      dataSource="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/173/515/GCA_002173515.1_ASM217351v1/GCA_002173515.1_ASM217351v1_genomic.gff.gz")

# Read your GFF file line by line
gff_lines <- readLines("./Genome/Bacterial_rna_seq/GCA_002173515.1_ASM217351v1_genomic.gff")

# Initialize an empty vector to store the gene IDs
gene_ids <- character(0)

# Define a regular expression pattern to match lines containing gene IDs
pattern <- "ID=gene-S101446_\\d+;"

# Loop through the lines and extract gene IDs
for (line in gff_lines) {
  # Check if the line contains the pattern
  if (grepl(pattern, line)) {
    # Extract the gene ID using regular expressions
    gene_id <- regmatches(line, regexec(pattern, line))
    if (!is.na(gene_id[[1]])) {
      # Remove "ID=gene-" and the trailing semicolon to get the gene ID
      gene_id <- gsub("ID=gene-|;", "", gene_id[[1]])
      gene_ids <- c(gene_ids, gene_id)
    }
  }
}

# Print the extracted gene IDs
cat(gene_ids, sep = "\n")



# Get the exons per gene, and compute bp lengths of all genes
exons <- exonsBy(db2, by = "gene")

metadata(exons)$gene_ids <- gene_ids 
head(metadata(exons)$gene_ids)
gene.lengths <- lapply(exons, function(x){ sum(width(reduce(x))) })






# List your BAM files
bam_files <- c(
  "./Genome/Bacterial_rna_seq/Data/Output/SRR20958305.aln/SRR20958305.sorted.RD.RG.bam",
  "./Genome/Bacterial_rna_seq/Data/Output/SRR20958306.aln/SRR20958306.sorted.RD.RG.bam",
  "./Genome/Bacterial_rna_seq/Data/Output/SRR20958307.aln/SRR20958307.sorted.RD.RG.bam",
  "./Genome/Bacterial_rna_seq/Data/Output/SRR20958308.aln/SRR20958308.sorted.RD.RG.bam",
  "./Genome/Bacterial_rna_seq/Data/Output/SRR20958309.aln/SRR20958309.sorted.RD.RG.bam",
  "./Genome/Bacterial_rna_seq/Data/Output/SRR20958310.aln/SRR20958310.sorted.RD.RG.bam"
  )

# Load BAM files
bfl <- BamFileList(bam_files)



# Perform read count summarization
se <- summarizeOverlaps(
  features = exons,
  reads = bfl,
  mode = "Union",
  singleEnd = FALSE,
  ignore.strand = TRUE,
  fragments = TRUE
)


# Extract the raw-reads per gene
readcount <- assay(se)
colnames(readcount) <- sub(".bam", "", colnames(readcount))

# Optionally, save the read counts to a file
write.table(readcount, "readcount.txt", sep = "\t", quote = FALSE)

# Get the total number of reads per sample
total_reads <- colSums(readcount)
print(readcount)
# Calculate FPKM (Fragments Per Kilobase per Million mapped fragments)
fragment_lengths <- rep(mean(width(exons) / 1000), length(total_reads))
fpkm <- (10^9 * readcount) / (total_reads * fragment_lengths)

# save the FPKM values to a file
write.table(fpkm, "fpkm.txt", sep = "\t", quote = FALSE)





# Normalize the data
normalize_data <- function(data) {
  normalized <- round(normalize.quantiles(as.matrix(data)), d = 1)
  colnames(normalized) <- colnames(data)
  rownames(normalized) <- rownames(data)
  return(normalized)
}





# Perform data normalization
fpkm.norm <- normalize_data(fpkm)




# save the RPKM values to a file
write.table(fpkm.norm, "fpkm.txt", sep = "\t", quote = FALSE)


 

# P-values and Log2 fold change
pvals <- apply(fpkm, 1, function(x) {
  # Use your statistical test for differential expression here
  tryCatch(t.test(x[1:3], x[4:6])$p.value, error = function(x) { return(NA); })
})

fc <- apply(fpkm, 1, function(x) {
  # Calculate the fold change based on your experimental design
  tryCatch(log2(mean(x[1:3]) / mean(x[4:6])), error = function(x){return(NA);})
})

# Assign colors based on P-values
colz <- rep("black", length(pvals))
colz[which(pvals < 5e-2)] <- "red"
colz[which(pvals < 1e-2)] <- "gold"
colz[which(pvals < 1e-3)] <- "blue"

# Volcano plot (x = fc, y = -log10(P-values))
plot(fc, -log10(pvals), col = colz, pch = 18, main = "Volcano Plot", xlab = "Fold Change")
legend("topleft", pch = 18, c("<0.05", "<0.01", "<0.001"), col = c("red", "gold", "blue"))

fpkm[is.na(fpkm)] <- 0
# Sample names from colnames of fpkm
sample_names <- colnames(fpkm)

# Create a DataFrame for colData
colData <- DataFrame(
  sampleName = sample_names,
  condition = factor(c("mid","mid", "late","late", "initial","initial")),
  replicate = factor(c(2,2,2))  # Example of replicate information
)
# If you have FPKM values, you can convert them to counts
counts <- round(fpkm * total_reads * fragment_lengths / 10^9)
counts[is.na(counts)] <- 0

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- DESeq(dds)
design_formula <- ~condition

dds <- DESeq(dds, test = "Wald", fitType = "parametric", minReplicatesForReplace = Inf)
res_mid_initial <- results(dds, contrast = c("condition", "mid", "initial"))
res_late_initial <- results(dds, contrast = c("condition", "late", "initial"))
res_mid_late <- results(dds, contrast = c("condition", "mid", "late"))

#Visualization of volcano plot

create_volcano_plot <- function(results, padj_threshold, log2fc_threshold, comparison_name) {
  # Create a data frame from DESeq2 results
  results_df <- as.data.frame(results)
  
  # Filter results based on thresholds
  significant_genes <- with(results_df, padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)
  
  # Create a volcano plot
  volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant_genes)) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "blue") +
    labs(title = paste("Volcano Plot:", comparison_name)) +
    theme_minimal()
  
  # Plot the volcano plot
  print(volcano_plot)
}

# Example usage:
create_volcano_plot(res_mid_late, padj_threshold = 0.05, log2fc_threshold = 1, comparison_name = "mid vs. late")
create_volcano_plot(res_mid_initial, padj_threshold = 0.05, log2fc_threshold = 1, comparison_name = "mid vs. initial")
create_volcano_plot(res_late_initial, padj_threshold = 0.05, log2fc_threshold = 1, comparison_name = "late vs. initial")



# Define a function for plotting upregulated 
plot_differentially_expressed_genesup <- function(results, n_genes, title, color) {
  # Extract gene names and log2 fold changes
  top_genes <- head(results[order(-results$log2FoldChange), ], n = n_genes)
  
  # Create a data frame for plotting
  top_genes_df <- data.frame(
    GeneName = rownames(top_genes),
    log2FoldChange = top_genes$log2FoldChange
  )
  
  
  
  # Plot the genes
  ggplot(top_genes_df, aes(x = GeneName, y = log2FoldChange)) +
    geom_bar(stat = "identity", fill = color) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
}
# Define a function for plotting downregulated genes
plot_differentially_expressed_genesdown <- function(results, n_genes, title, color) {
  # Extract gene names and log2 fold changes
  top_genes <- head(results[order(results$log2FoldChange), ], n = n_genes)
  
  # Create a data frame for plotting
  top_genes_df <- data.frame(
    GeneName = rownames(top_genes),
    log2FoldChange = top_genes$log2FoldChange
  )
  
  # Plot the genes
  ggplot(top_genes_df, aes(x = GeneName, y = log2FoldChange)) +
    geom_bar(stat = "identity", fill = color) +
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
}
# res_mid_initial upregulated
plot_differentially_expressed_genesup(res_mid_initial, n_genes = 10, title = "Top Upregulated Genes mid_initial", color = "red")
# res_mid_initial downregulated
plot_differentially_expressed_genesdown(res_mid_initial, n_genes = 10, title = "Top Downregulated Genes mid_initial", color = "blue")

# res_mid_late upregulated
plot_differentially_expressed_genesup(res_mid_late, n_genes = 10, title = "Top Upregulated Genes mid_late", color = "orange")
# res_mid_initial downregulated
plot_differentially_expressed_genesdown(res_mid_late, n_genes = 10, title = "Top Downregulated Genes mid_late", color = "purple")

# res_late_initial upregulated
plot_differentially_expressed_genesup(res_late_initial, n_genes = 10, title = "Top Upregulated Genes late_initial", color = "darkgreen")
# res_late_initial downregulated
plot_differentially_expressed_genesdown(res_late_initial, n_genes = 10, title = "Top Downregulated Genes late_initial", color = "#FF5733")
