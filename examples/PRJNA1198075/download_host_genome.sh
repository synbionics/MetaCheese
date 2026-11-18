#!/usr/bin/env bash
# Download the host genome of the examples, Bos taurus ARS-UCD2.0
# (GCA_002263795.4), from NCBI into data/gene/Bos_taurus/, and check it.
# The pipeline builds the Bowtie2 index from it on the first run (about 40 min
# with 16 threads), or build it now with:
#   docker run --rm -v "$PWD/data:/main/data" metacheese \
#       bash scripts/build_bowtie2_index.sh Bos_taurus 16
#
# Usage, from the repository root: bash examples/PRJNA1198075/download_host_genome.sh
set -euo pipefail

ACC=GCA_002263795.4_ARS-UCD2.0
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/263/795/$ACC
DIR=data/gene/Bos_taurus
FASTA=$DIR/${ACC}_genomic.fna
# MD5 of the uncompressed FASTA the tested index was built from (2026-11)
FASTA_MD5=01a077ef1e8cc0e961bc995fe1791319

[[ -d data/gene ]] || { echo "ERROR: run from the repository root (data/gene/ not found)" >&2; exit 1; }
mkdir -p "$DIR"

if [[ -s "$FASTA" ]] && echo "$FASTA_MD5  $FASTA" | md5sum -c --quiet 2>/dev/null; then
    echo "[OK] $FASTA already present and checked"
    exit 0
fi

tmp="$DIR/${ACC}_genomic.fna.gz.part"
echo "Downloading $ACC (about 850 MB) ..."
curl -sSfL --retry 3 -o "$tmp" "$URL/${ACC}_genomic.fna.gz"

# NCBI's checksum of the compressed file, then ours of the FASTA
gz_md5=$(curl -sSfL "$URL/md5checksums.txt" | tr -d '\r' | awk -v f="./${ACC}_genomic.fna.gz" '$2 == f {print $1}')
[[ -n "$gz_md5" ]] || { echo "ERROR: no checksum for ${ACC}_genomic.fna.gz in md5checksums.txt" >&2; exit 1; }
echo "$gz_md5  $tmp" | md5sum -c --quiet || { echo "ERROR: download corrupted, run again" >&2; rm -f "$tmp"; exit 1; }

gunzip -c "$tmp" > "$FASTA.part"
rm -f "$tmp"
echo "$FASTA_MD5  $FASTA.part" | md5sum -c --quiet \
    || { echo "ERROR: the FASTA differs from the one the pipeline was tested with" >&2; exit 1; }
mv "$FASTA.part" "$FASTA"
echo "[OK] $FASTA"
