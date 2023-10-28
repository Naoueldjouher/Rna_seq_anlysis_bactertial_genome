#Fetch the reference genome (replace with your organism and ID)
./edirect/efetch -db nucleotide -id CP021467.1 -format fasta
#Download the reference GTF for your organism and compress it
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/173/515/GCA_002173515.1_ASM217351v1/GCA_002173515.1_ASM217351v1_genomic.gtf.gz
gunzip -c GCA_002173515.1_ASM217351v1_genomic.gtf.gz > GCA_002173515.1_ASM217351v1_genomic.gtf

# Index the reference genome using samtools (replace with your organism and file)
samtools faidx CP021467.1.fa

# Generate genome/transcriptome index using Bowtie
 Bowtie2-build  -i ./CP021467.1.fa.gz  -o ./Bowtie2/CP021467.1_index
