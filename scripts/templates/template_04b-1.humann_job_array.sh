#!/bin/bash
# Step 04b-1 (diagram C.3): HUMAnN functional profiling, one sample at a time.
set -euo pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate humannenv

NUC_DB="/main/db/humann/chocophlan/chocophlan"
PROT_DB="/main/db/humann/uniref/uniref"

# Check che i db ci siano
if [ ! -d "$NUC_DB" ] || [ -z "$(ls -A "$NUC_DB")" ]; then
    echo "ERRORE: DB ChocoPhlAn mancante" >&2; exit 1
fi
if [ ! -d "$PROT_DB" ] || [ -z "$(ls -A "$PROT_DB")" ]; then
    echo "ERRORE: DB UniRef90 mancante" >&2; exit 1
fi


# Percorsi input/output (lascia come i tuoi)
taxon_dir="@04b_1_var1@"
input_dir="@04b_1_var2@"
output_dir="@04b_1_var3@"

mkdir -p "$output_dir"

# HUMAnN uses the taxonomic profiles from step 04 (MetaPhlAn 4.1.0, bioenv), but
# still checks the MetaPhlAn version before starting. The metaphlan 4.2.x that
# conda installs in humannenv prints an extra "No complete MetaPhlAn Bowtie2
# database found" line, which breaks that check ("Can not call software version
# for metaphlan"). Point HUMAnN at the bioenv MetaPhlAn, the one that produced
# the profiles, via a folder holding only that executable (adding all of
# bioenv/bin to PATH would also shadow HUMAnN's bowtie2/diamond).
MPA_BIN="$(mktemp -d)"
ln -s /opt/conda/envs/bioenv/bin/metaphlan "$MPA_BIN/metaphlan"
trap 'rm -rf "$MPA_BIN"' EXIT

# Loop sui campioni
samples=($(ls -1 ${input_dir}/*fq.1.gz | sort))
for sample in "${samples[@]}"; do
    nome_campione=$(basename "$sample" .fq.1.gz)
    input_file="${input_dir}/${nome_campione}.fq.1.gz"
    taxonomic_profile="${taxon_dir}/${nome_campione}.txt"

    # Resume support: a sample whose three output tables exist is already done
    # (HUMAnN writes them only at the end; an interrupted sample leaves just a
    # *_humann_temp_* folder, which is removed and recomputed).
    done_prefix="$output_dir/${nome_campione}.fq"
    if [[ -s "${done_prefix}_genefamilies.tsv" && -s "${done_prefix}_pathabundance.tsv" && -s "${done_prefix}_pathcoverage.tsv" ]]; then
        echo "Sample $nome_campione already done, skipping"
        continue
    fi
    rm -rf "${done_prefix}"_humann_temp_*

    if [[ -f "$input_file" && -f "$taxonomic_profile" ]]; then
        echo "Processing sample: $nome_campione"
        humann \
          --input "$input_file" \
          --output "$output_dir" \
          --threads @04b_1_par1@ \
          --memory-use maximum \
          --diamond-options "--top 1 --outfmt 6 --block-size @04b_1_par2@" \
          --remove-temp-output \
          --taxonomic-profile "$taxonomic_profile" \
          --metaphlan "$MPA_BIN" \
          --nucleotide-database "$NUC_DB" \
          --protein-database "$PROT_DB"
    else
        echo "WARNING: Input or taxonomic profile missing for sample: $nome_campione. Skipping..." >&2
    fi

    echo "HUMAnN processing complete for sample: $nome_campione"
done
