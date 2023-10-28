#!/bin/bash


input_base="$1"


for input_base in "${input_bases[@]}"; do
    index_path="./Genome/Bacterial_rna_seq/Bowtie2/CP021467.1_index"
    output_dir="./Genome/Bacterial_rna_seq/Data/Output/$input_base.aln"
    bowtie2_exec="/usr/local/bin/bowtie2"
    samtools="/opt/anaconda3/envs/myenv2/bin/samtools"
    echo "Processing input base: $input_base"

    # Run bowtie2
    "$bowtie2_exec" -x "$index_path" \
        -1 "$output_dir/$input_base"_1.P.fastq.gz \
        -2 "$output_dir/$input_base"_2.P.fastq.gz \
        -S "$output_dir/$input_base.sam"

    # Convert SAM to BAM
    "$samtools" view -bS "$output_dir/$input_base.sam" \
        -o "$output_dir/$input_base.bam"

    # Sort the BAM file
    "$samtools" sort "$output_dir/$input_base.bam" \
        -o "$output_dir/$input_base.sorted.bam"

    # Index the sorted BAM file
    "$samtools" index "$output_dir/$input_base.sorted.bam"

    # Run samtools flagstats and coverage
    "$samtools" flagstat "$output_dir/$input_base.sorted.bam"
    "$samtools" coverage "$output_dir/$input_base.sorted.bam"
done
