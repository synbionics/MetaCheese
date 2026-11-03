#!/usr/bin/env bash
set -euo pipefail

# Activate Conda (you can also use conda run directly in commands)
source /opt/conda/etc/profile.d/conda.sh
conda activate bioenv

# Metaphlan database path (absolute and consistent with the Dockerfile)
DEFAULT_DB_FOLDER="/main/db/metaphlan"

# Sample and output folders (replaced by Nextflow)
sample_folder="@04_var2@"
output_folder="@04_var3@"
mkdir -p "$output_folder"

# Check that the database exists
if [ ! -d "$DEFAULT_DB_FOLDER" ]; then
    echo "Error: DB folder not found → $DEFAULT_DB_FOLDER" >&2
    exit 1
fi

# Loop through all samples
for file1 in "$sample_folder"/*.fq.1.gz; do
    nome_campione=$(basename "$file1" .fq.1.gz)
    file2="$sample_folder/${nome_campione}.fq.2.gz"

    if [[ ! -f "$file2" ]]; then
        echo "WARNING: sample $nome_campione - file .fq.2.gz missing, skipping."
        continue
    fi

    echo "Running Metaphlan analysis for sample $nome_campione..."

    metaphlan \
        "$file1","$file2" \
        --input_type fastq \
        --bowtie2db "$DEFAULT_DB_FOLDER" \
        --index mpa_vJun23_CHOCOPhlAnSGB_202307 \
        --nproc @04_par1@ \
        --bowtie2out "$output_folder/${nome_campione}.bowtie2.bz2" \
        -s "$output_folder/${nome_campione}.sam.bz2" \
        -o "$output_folder/${nome_campione}.txt"
done

# Merge taxonomic profiles if results are available
txt_count=$(ls "$output_folder"/*.txt 2>/dev/null | wc -l)
if [[ "$txt_count" -eq 0 ]]; then
    echo "No .txt files found in $output_folder: skipping merge_metaphlan_tables.py"
else
    echo "Merging taxonomic profiles..."
    merge_metaphlan_tables.py "$output_folder"/*.txt > "$output_folder/merged_abundance_table.txt"
fi

# Alpha diversity analysis
for metric in richness shannon simpson gini; do
    echo "Calculating alpha diversity ($metric)..."
    Rscript /main/data/calculate_diversity.R \
        -f "$output_folder/merged_abundance_table.txt" \
        -d alpha \
        -m "$metric" \
        -o "$output_folder/diversity"
done

# Beta diversity analysis (aitchison)
for metric in bray-curtis jaccard clr ; do
    echo "Calculating beta diversity ($metric)..."
    Rscript /main/data/calculate_diversity.R \
        -f "$output_folder/merged_abundance_table.txt" \
        -d beta \
        -m "$metric" \
        -o "$output_folder/diversity"
done

TREE_PATH="$output_folder/tree.nwk"
if [[ ! -f "$TREE_PATH" ]]; then
    echo "WARNING: phylogenetic tree tree.nwk missing → skipping UniFrac metrics!"
fi

# Run UniFrac metrics only if the tree is present
if [[ -f "$TREE_PATH" ]]; then
    for metric in weighted-unifrac unweighted-unifrac; do
        echo "Calculating beta diversity ($metric)..."
        Rscript /main/data/calculate_diversity.R \
            -f "$output_folder/merged_abundance_table.txt" \
            -d beta \
            -m "$metric" \
            -t "$TREE_PATH" \
            -o "$output_folder/diversity"
    done
else
    echo "Skipping UniFrac metrics (weighted/unweighted) due to missing tree.nwk."
fi

conda deactivate
