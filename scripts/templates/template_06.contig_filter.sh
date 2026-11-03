#!/bin/bash
set -euo pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate bioenv

cd @06_var1@
mkdir -p @06_var2@

echo "Copying contigs to @06_var2@"
find . -type f -name "contigs.fasta" -print0 | while IFS= read -r -d '' file; do
    parent_dir_name=$(basename "$(dirname "$file")")
    cp "$file" "@06_var2@/${parent_dir_name}.fasta"
done

cd @06_var2@
mkdir -p @06_var3@

echo "Filtering contigs in @06_var3@ keeping sequences >500 bp"

find -type f -name "*.fasta" -print0 | while IFS= read -r -d '' fasta_file; do
    sample_name=$(basename "$fasta_file" .fasta)
    # Measure the sequence only: with RS=">" each record also holds the header
    # line and the line breaks, so length($0) would let shorter contigs through.
    awk 'BEGIN{RS=">"; ORS=""} NR > 1 {
        i = index($0, "\n"); seq = substr($0, i + 1); gsub(/\n/, "", seq)
        if (length(seq) > 500) print ">"$0
    }' "$fasta_file" > "@06_var3@/${sample_name}.fasta_sort.fasta"
done

echo "Filtering completed."
conda deactivate
