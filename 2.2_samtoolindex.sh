#!/bin/bash

# Define an array of input bases (assuming space-separated)
input_bases=($1)

# Path to the samtools executable
samtools="/opt/anaconda3/envs/myenv2/bin/samtools"

# Java command for Picard
java_cmd="java -Xmx4g -jar"

# Path to the Picard tools executable
picard="/opt/anaconda3/envs/myenv2/share/picard-2.27.5-0/picard.jar"

# Define the base output directory
output_base_dir="./Genome/Bacterial_rna_seq/Data/Output"

for input_base in "${input_bases[@]}"; do
    output_dir="${output_base_dir}/${input_base}.aln"
    outbase="${output_dir}/${input_base}"
    bam_file="${outbase}.sorted.bam"
    p_bam="${outbase}.sorted.RD.bam"

    # STEP 3.1 - Create a samtools index
    $samtools index "$p_bam"

    # STEP 3.2 - Create mapping and coverage statistics
    $samtools flagstat "$p_bam"
    $samtools coverage "$p_bam"

    # STEP 4 - Add read group (1) and sample run, library, and name
    rg_bam="${output_dir}/${input_base}.sorted.RD.RG.bam"
    rg_opts="-PL ILLUMINA -PU run -LB $(echo $input_base | sed 's/SRR//') -SM $input_base"
    p_cmd="$java_cmd $picard AddOrReplaceReadGroups -I $p_bam -O $rg_bam $rg_opts"
    $p_cmd

    # STEP 4.1 - Create a samtools index for the final BAM file
    $samtools index "$rg_bam"
done
