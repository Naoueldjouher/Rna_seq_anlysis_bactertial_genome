#Set up 
#Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
#Conda 
brew install --cask miniconda
#R
conda create -n myenv r-base
#wget
brew install wget
#Java
brew tap AdoptOpenJDK/openjdk
brew cask install adoptopenjdk11
#edirect
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
#Add edirect to Your Shell Environment
export PATH=$PATH:$HOME/edirect
#Replace $HOME/edirect with the actual path where edirect is installed.
#Load the updated configuration with:
source ~/.bashrc

# Create a directory for software and navigate to it
mkdir software
cd software

# Clone Trimmomatic repository and build it
git clone https://github.com/usadellab/Trimmomatic.git

# Download and extract SRA Toolkit
cd software
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-centos_linux64-cloud.tar.gz
mkdir sratoolkit
tar -xzf sratoolkit.3.0.0-centos_linux64-cloud.tar.gz –C sratoolkit
./sratoolkit/usr/local/ncbi/sra-tools/bin/vdb-config --interactive
#Install htslib,samtools,bcftools using conda
conda install -c bioconda htslib
conda install -c bioconda samtools
conda install -c bioconda bcftools
#Install Bowtie2 using HomeBrew
cd software
brew install bowtie2
# Download and extract IGV (Integrative Genomics Viewer)
cd software
wget https://data.broadinstitute.org/igv/projects/downloads/2.14/IGV_Linux_2.14.1_WithJava.zip
unzip IGV_Linux_2.14.1_WithJava.zip
#Install multiqc,fastqc
conda install -c bioconda multiqc
conda install -c bioconda fastqc
#In R 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicFeatures")
BiocManager::install("preprocessCore")
BiocManager::install("RColorBrewer")
BiocManager::install("rtracklayer")
BiocManager::install("BiocParallel")
BiocManager::install("GenomicAlignments")
BiocManager::install("DESeq2")


