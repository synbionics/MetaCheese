#!/bin/bash

# Activate conda
source /opt/conda/etc/profile.d/conda.sh
conda activate bowtieenv

BASE="data/gene"

# Usage: bash scripts/build_bowtie2_index.sh [SPECIES] [THREADS]
# Species folder name (e.g. Bos_taurus): from the first argument, or asked
SPECIES="${1:-}"
THREADS="${2:-4}"
[[ -n "$SPECIES" ]] || read -p "Species folder name inside $BASE/: " SPECIES
DIR="$BASE/$SPECIES"

# Check if folder exists
if [ ! -d "$DIR" ]; then
  echo "[ERROR] Folder $DIR does not exist."
  exit 1
fi

# Check for FASTA files
FASTA=$(ls "$DIR"/*.fa "$DIR"/*.fna "$DIR"/*.fasta 2>/dev/null | head -n1)
if [ -z "$FASTA" ]; then
  echo "[ERROR] No FASTA file (.fa/.fna/.fasta) found in $DIR."
  exit 1
fi

# Index prefix (same name as FASTA without extension)
PREFIX="$DIR/$SPECIES"

# Check if index already exists
if ls "${PREFIX}"*.bt2 &>/dev/null; then
  # non-interactive call (species given as argument): keep the existing index
  [[ -n "${1:-}" ]] && { echo "[OK] Index already exists: $PREFIX"; conda deactivate; exit 0; }
  read -p "Index already exists, do you want to rebuild it? [y/N]: " ans
  [[ "$ans" != "s" && "$ans" != "S" && "$ans" != "y" && "$ans" != "Y" ]] && { echo "Exiting."; conda deactivate; exit 0; }
fi

# Build index
echo "[INFO] Building index for $FASTA ..."
bowtie2-build --threads "$THREADS" "$FASTA" "$PREFIX"

conda deactivate
echo "[OK] Index created with prefix: $PREFIX"
