#!/bin/bash
set -euo pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate spadesenv

campioni_folder="@05_var1@"
output_folder="@05_var2@"
mkdir -p "$output_folder"

for fq1 in "$campioni_folder"/*.fq.1.gz; do
    nome_campione=$(basename "$fq1" .fq.1.gz)
    fq2="$campioni_folder/${nome_campione}.fq.2.gz"

    if [[ ! -f "$fq2" ]]; then
        echo "WARNING: File $fq2 missing for $nome_campione, it will be skipped." >&2
        conda deactivate
        exit 1
    fi

    echo "Processing sample: $nome_campione"

    # SPAdes 4 accepts only known read extensions (.fq.gz, .fastq.gz, ...),
    # not step 03's <sample>.fq.1.gz, so it reads them through links.
    reads_dir="$output_folder/reads_links/$nome_campione"
    mkdir -p "$reads_dir"
    ln -sfn "$(realpath "$fq1")" "$reads_dir/R1.fq.gz"
    ln -sfn "$(realpath "$fq2")" "$reads_dir/R2.fq.gz"

    if ! spades.py -1 "$reads_dir/R1.fq.gz" -2 "$reads_dir/R2.fq.gz" \
        --meta -t @05_par1@ --memory @05_par2@ \
        --only-assembler \
        -o "$output_folder/$nome_campione"; then
        echo "ERROR: SPAdes failed for $nome_campione" >&2
        exit 1
    fi

    echo "Sample completed: $nome_campione"
done

echo "All valid samples completed."
conda deactivate
