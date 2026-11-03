#!/bin/bash
set -euo pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate bowtieenv

output_folder="@07_var1@"
input_folder="@07_var2@"
mkdir -p "$output_folder"

echo "Creating Bowtie2 indexes in $output_folder from fasta files in $input_folder"

found=false
for file in "$input_folder"/*.fasta; do
    [[ -e "$file" ]] || continue
    found=true

    base_name=$(basename "$file" .fasta)
    echo "Indexing $base_name"

    if ! bowtie2-build --threads @07_par1@ "$file" "$output_folder/${base_name}_index_base"; then
        echo "ERROR: bowtie2-build failed for $file" >&2
        exit 1
    fi
done

if ! $found; then
    echo "No .fasta file found in $input_folder" >&2
    exit 1
fi

echo "Index creation completed successfully."
conda deactivate
