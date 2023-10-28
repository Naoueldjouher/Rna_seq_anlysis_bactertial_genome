<h1>RNA-Seq Workflow Analysis for Komagataeibacter Europaeus</h1>

<p>
 This study utilizes RNA-Seq analysis to investigate the mechanisms of acid resistance (AR) in Acetic Acid Bacteria (AAB) during industrial vinegar fermentation. It explores the influence of various fermentation stages (initial, mid, and final) and potential physiological adaptations on AR.





</p>

<p>
  <strong>Research Question:</strong><br>
Which genes display significant upregulation and downregulation in response to increased acidity during industrial vinegar fermentation?


<p>
  <strong>Link to the Study:</strong><br>
  <a href="https://www.frontiersin.org/articles/10.3389/fmicb.2022.956729/ful" target="_blank">Read the full study</a>
</p>

<p>
  <strong>Inspired by:</strong><br>
  Dr. Arends' RNA-Seq analysis pipeline (<a href="https://www.youtube.com/watch?v=PlqDQBl22DI&list=PLhR2Go-lh6X63hnyBzwWNvsaw1R79ESPI&pp=iAQB" target="_blank">Watch here</a>
</p>
<p>
<h2>Workflow:</h2>
<div style="text-align: center;">
  <img src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/cc112a42-8765-4988-a367-a3f01f7c5555" alt="workflow rna" width="400">
</div>

 
</p>

  <p><strong>Step 1:</strong> Install the necessary software for your RNA-Seq analysis.</p>
  <p><strong>Step 2:</strong> Obtain the FASTA file of the reference genome and the GTF files. Index the reference genome for efficient analysis.</p>
  <p><strong>Step 3:</strong> Download the SRR files containing your sequencing data. Check the quality of the data using FastQC.</p>
  <p><strong>Step 4:</strong> Trim and preprocess your sequencing data using Trimmomatic for quality control.</p>
  <p><strong>Step 5:</strong> Perform sequence alignment to map the reads to the reference genome.</p>
  <p><strong>Step 6:</strong> Utilize Picard for  manipulations and analyses of BAM files, particularly for read grouping .</p>
  <p><strong>Step 7:</strong> Index the processed data using Samtools for subsequent analysis and visualization.</p>
  <p><strong>Step 8:</strong> Perform normalization procedures to prepare your data for differential expression analysis.</p>
  <p><strong>Step 9:</strong> Extract differentially expressed genes and visualize the results to gain insights into your RNA-Seq experiment.</p>


<p>

<h2>Installation of Tools:</h2>


<h3>For Data Gathering, Cleaning, and Sequence Alignment:</h3>
<ol>
  <li><strong>Trimmomatic:</strong> For quality control and adapter trimming of raw sequencing data.</li>
  <li><strong>SRA Toolkit:</strong> To access and convert SRA data into FASTQ format.</li>
  <li><strong>Bowtie 2:</strong> A fast and efficient aligner used to map sequenced reads to a reference genome.</li>
  <li><strong>Picard:</strong> Essential for manipulating and analyzing BAM files, particularly for read grouping and sequence dictionary creation.</li>
  <li><strong>Samtools:</strong> Provides utilities for working with sequence alignment data in the SAM/BAM format.</li>
  <li><strong>Htslib:</strong> The foundational library used by Samtools and Bcftools for accessing high-throughput sequencing data.</li>
  <li><strong>Bcftools:</strong> Primarily used for variant calling and manipulation of VCF/BCF files.</li>
  <li><strong>IGV (Integrative Genomics Viewer):</strong> A visualization tool for exploring and interpreting genomic data.</li>
  <li><strong>FastQC:</strong> FastQC is used for quality control of raw sequencing data.</li>
  <li><strong>MultiQC:</strong> for aggregating and summarizing results.</li>
</ol>
<h3>For DEG analysis and Data visualization:</h3>
<ol>
  <li><strong>Rsamtools:</strong> Provides interface to the SAM/BAM file format.</li>
  <li><strong>GenomicFeatures:</strong> Representation of genomic annotations.</li>
  <li><strong>PreprocessCore:</strong> Offers functions for data preprocessing and quality control in high-throughput sequencing data.</li>
  <li><strong>RColorBrewer:</strong> Provides additional color palettes for creating visually appealing plots and visualizations.</li>
  <li><strong>rtracklayer:</strong> Allows reading and writing data in common NGS file formats.</li>
  <li><strong>BiocParallel:</strong> Enables parallel processing and efficient computation in Bioconductor packages.</li>
  <li><strong>GenomicAlignments:</strong> Handles the alignment of high-throughput sequencing data to the reference genome.</li>
  <li><strong>DESeq2:</strong> Performs differential gene expression analysis of RNA-seq data.</li>
  <li><strong>ggplot:</strong> Data visualization.</li>
</ol>
<h2>Data Gathering, Cleaning, and Alignment:</h3>

<h3>Setting up the required data for RNA-seq experiment:</h3>

<p> - There was no vcf file available</p>
<p> - Creation and Building an index of the reference genome using Bowtie2 enables the quick alignment of read to the reference genome Here, you are using Bowtie2 to build an index for your reference genome.</p>
<p> The reference genome can be downloaded from here: <a href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/173/515/GCA_002173515.1_ASM217351v1/">Download Link</a> </p>
<p> The SRR numbers can be downloaded from here: <a href="https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP390398&o=acc_s%3Aa">Download Link</a> </p>
<p>-There are 6 runs each 2 were observed at a different stage (final, mid, initial)</p>
<p><strong>-<span style="color:#FF5733;">SRR20958305</span></strong> and <strong><span style="color:#FF5733;">SRR20958306</span></strong>  during the final stage</p>
<p><strong>-<span style="color:#FF5733;">SRR20958307</span></strong> and <strong><span style="color:#FF5733;">SRR20958308</span></strong>  during the mid stage</p>
<p><strong>-<span style="color:#FF5733;">SRR20958309</span></strong> and <strong><span style="color:#FF5733;">SRR20958310</span></strong> during the initial stage</p>
  <h3>1. Fastq check using Multiqc:</h3>

<div style="display: flex; justify-content: center;">
  <img width="627" alt="Screenshot 2023-10-28 at 20:52:28" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/ebf98d37-1802-4ee5-a92d-829481944110">
  <img width="627" alt="Screenshot 2023-10-28 at 20:52:28" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/e683c915-9b89-433f-9d46-d4afa1210bf4">
  <img width="627" alt="Screenshot 2023-10-28 at 20:52:28" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/7537861d-fc9e-443b-a61e-f56e62ba9f7c">
</div>
<p> -The overall data quality score is good.</p>
<p> - A significant proportion of duplicated reads is noticeable, this can often result from various sources, including PCR amplification biases or happens when dealing with bacterial genome  these duplicated reads may not contribute additional biological information and can artificially inflate the apparent significance of certain features .
</p>
<h3>Visualization of results from Picard:</h3>
<img width="627" alt="Screenshot 2023-10-28 at 21 56 52" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/596484ce-6be7-413e-a3ae-1c0270b9733b">

<p>Trimming and the use of Picard tools have been effective in reducing the percentage of duplicated paired reads, although some level of duplication persists. 
It's worth noting that the presence of unmapped regions in the data could be attributed to the nature of bacterial genomes. Bacterial genomes often contain repetitive elements, mobile genetic elements, and genomic rearrangements,Additionally, bacterial genomes may exhibit variations in regions that are not yet fully characterized or annotated. 
To address this,downstream analyses, such as variant calling and annotation, should account for potential genomic variations and unmapped regions,which was not possible to perform for this study.
</p>
<p>Overrepresented regions in RNA-Seq data of a bacterial genome can arise from various sources, including sequencing artifacts, biological variability, and experimental techniques like PCR amplification.</p>

<h3>Visualization of the diffentially expressed genes</h3>
<h4>Volcano plot of the differentially expressed genes during the late,mid,initial</h4>
<div style="display: flex; justify-content: space-between;">
  <img width="250" alt="Rplot01" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/8451b39e-8f1b-4476-9b13-96f2bf8f1565">
  <img width="250" alt="Rplot03" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/f7c2c567-e06a-4cd3-b29e-0e3711ef711b">
  <img width="226" alt="mid_initial" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/8bee15d0-7a75-4da8-a1d1-dba64de8ed5b">
</div>
<p><strong>1-mid vs late</strong></p>
<p>In the "mid vs. late" comparison, only a limited number of genes exhibit significance. Among these, a majority are upregulated, indicating that the gene expression changes during the transition from the "mid" to the "late" condition are primarily characterized by an increase in expression levels.</p>
<p><strong>2-mid vs. Initial:</strong></p>
<p>In the "mid" vs. "initial" comparison, you observe a higher number of significant genes. These genes are scattered throughout the plot, and there's a mix of upregulated and downregulated genes. This suggests a more complex and dynamic gene expression response in the "mid" condition compared to the "initial."</p>
<p><strong>3-Late vs. Initial:</strong></p>
<p>In the "late" vs. "initial" comparison, you see relatively few significant genes, and they appear to be aligned in a line. This alignment suggests that these genes might be part of a specific pathway or process that is triggered during the transition from the "initial" to the "late" condition. The linear arrangement of significant genes indicates a coordinated response.</p>
<h4>Differentially upregulatted genes</h4>
<div style="display: flex; justify-content: space-between;">
  <img width="250" alt="upregulated mid_initial" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/ef3c47b0-c493-4d0d-9ee3-2d9d207d02a3">
  <img width="250" alt="upregulated mid_Late" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/9461729f-7152-4a4e-8d0c-b1f7bba80ee2">
  <img width="226" alt="top_regulated late_initial" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/546e0bfb-b441-4cd2-b550-231b9e1d00d1">
</div>
<p>-For mid vs initial genes <strong>S101446_01327 and S101446_01326</strong>, including the vanA gene, exhibited significant upregulation. This pattern was consistent with the upregulation of <strong>101446_01655</strong> during the "mid_late" stage. Furthermore, these genes were found to be repeatedly upregulated in the "mid_initial" and "mid_late" transitions, suggesting their central role in adapting to these conditions.
</p>
<p>-Conversely,<strong> S101446_0063</strong> was upregulated in the "mid_late" transition and downregulated in "late_initial," indicative of a dynamic shift in gene expression. Similarly, <strong>S101446_00582</strong> displayed upregulation during the "late_initial" stage. These findings hint at the possible involvement of these genes in the adaptive responses to the respective conditions, which warrant further investigation.</p>
<h4>Differentially downregulatted genes</h4>
<div style="display: flex; justify-content: space-between;">
  <img width="250" alt="downregulated mid_initial" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/add1b8b6-6856-46f4-b84f-673b1cc73909">
  <img width="250" alt="downregulated mid_Late" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/1781736b-1346-4d7f-aaff-e90083d1f600">
  <img width="250" alt="down_regulated late_initial" src="https://github.com/Naoueldjouher/Rna_seq_anlysis_bactertial_genome/assets/80243706/db73395c-1a91-4aec-bbbd-63ab7ac11d7a">
</div>
<p>The highest expression peaks were observed in genes <strong>S101446_00682 and S101446_01866 </strong> for the case of late vs initial.In contrast, genes <strong>S101446_02564 and S101446_02771</strong> exhibited the highest expression levels during for the mid-initial . Additionally, genes <strong>S101446_02771 and S101446_02564</strong> were consistently downregulated when comparing "mid vs. initial" and "mid vs. late," indicating their consistent role in response to these conditions.These expression patterns suggest that certain genes play pivotal roles during specific transitions, highlighting their importance in adapting to different conditions within the experiment. </p>



