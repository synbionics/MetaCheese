#!/bin/bash
# Step 00 (diagram Q.1): quality check of the raw reads with FastQC.
# Produces one HTML report + one .zip (with fastqc_data.txt) per FASTQ file,
# then MultiQC merges them into a single page: multiqc_report.html.
# It only reports: reads are not modified, trimming happens in step 00-01.
set -euo pipefail
shopt -s nullglob

source /opt/conda/etc/profile.d/conda.sh
conda activate qcenv

input_dir="@00_dir1@"
output_dir="@00_dir2@"
threads="@00_par1@"
mkdir -p "$output_dir"

# Same naming conventions accepted by step 00-01: *_1/_2 or *_L1/_L2, .fq.gz
files=("$input_dir"/*_1.fq.gz "$input_dir"/*_2.fq.gz "$input_dir"/*_L1.fq.gz "$input_dir"/*_L2.fq.gz)
(( ${#files[@]} > 0 )) || { echo "ERROR: no *_1.fq.gz/*_2.fq.gz files in $input_dir" >&2; exit 1; }

echo "[$(date +'%F %T')] FastQC on ${#files[@]} files with $threads threads -> $output_dir"
# FastQC processes one file per thread
fastqc --threads "$threads" --outdir "$output_dir" --quiet "${files[@]}"

echo "[$(date +'%F %T')] FastQC completed"
conda deactivate

# One combined page for all samples (R1 and R2 of each sample side by side)
conda activate multiqcenv
multiqc --force --quiet --module fastqc \
    --title "MetaCheese raw-read quality (FastQC)" \
    --filename multiqc_report.html --outdir "$output_dir" "$output_dir"
conda deactivate
echo "[$(date +'%F %T')] MultiQC report: $output_dir/multiqc_report.html"
