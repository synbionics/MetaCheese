#!/bin/bash
set -euo pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate bowtieenv

FASTQ_DIR="@08-09_var1@"
INDEX_DIR="@08-09_var2@"
OUTPUT_DIR="@08-09_var3@"
mkdir -p "$OUTPUT_DIR"

for fq1 in "$FASTQ_DIR"/*.fq.1.gz; do
    base_name=$(basename "$fq1")
    file_name="${base_name%.fq.1.gz}"
    fq2="$FASTQ_DIR/${file_name}.fq.2.gz"

    if [[ ! -f "$fq2" ]]; then
        echo "WARNING: Mate file missing $fq2. Skipping $file_name." >&2
        exit 1
    fi

    index_base="${INDEX_DIR}/${file_name}.fasta_sort_index_base"
    if [[ ! -f "${index_base}.1.bt2" ]]; then
        echo "ERROR: Bowtie2 index missing for $file_name in $INDEX_DIR" >&2
        exit 1
    fi

    echo "Processing sample: $file_name"

    bam_file="${OUTPUT_DIR}/${file_name}.bam"
    sorted_bam="${OUTPUT_DIR}/${file_name}.sorted.bam"
    idxstats_file="${OUTPUT_DIR}/${file_name}.sorted.bam.idxstat"

    # Run Bowtie2 + samtools view
    bowtie2 -x "$index_base" \
        -1 "$fq1" -2 "$fq2" \
        -q --no-unal --very-sensitive-local \
        -p @08-09_par1@ 2> "${OUTPUT_DIR}/${file_name}_bowtie2.log" | \
        samtools view -bS -o "$bam_file" - 

    # Sort
    samtools sort -@ @08-09_par1@ "$bam_file" -o "$sorted_bam"

    # Index
    samtools index "$sorted_bam"

    # Stats
    samtools idxstats "$sorted_bam" > "$idxstats_file"

    echo "Sample $file_name completed"
done

echo "Mapping and statistics completed for all samples"
conda deactivate
