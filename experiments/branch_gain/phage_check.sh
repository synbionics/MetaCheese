#!/usr/bin/env bash
# What are a sample's high-copy contigs? Predicts proteins on the contigs above
# MIN_DEPTH and searches them against the image's UniRef90 with DIAMOND. Needs
# about 16 GB RAM. Then name_hits.py adds UniProt names (needs internet).
# This is evidence, not a verdict: run geNomad or CheckV before claiming phage.
#
# Usage: bash experiments/branch_gain/phage_check.sh <run_folder> <sample> <outdir>
# Environment: MIN_DEPTH (default 2000), THREADS (default 8), IMAGE (default metacheese)
set -euo pipefail
RUN="$(realpath "$1")"; S="$2"; OUT="$3"; mkdir -p "$OUT"; OUT="$(realpath "$OUT")"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MIN_DEPTH="${MIN_DEPTH:-2000}"; THREADS="${THREADS:-8}"; IMAGE="${IMAGE:-metacheese}"

awk -F'\t' -v d="$MIN_DEPTH" 'NR==FNR { if (FNR > 1 && $3 + 0 > d) hi[$1]; next }
    /^>/ { w = (substr($1, 2) in hi) } w' \
    "$RUN/10_metabat_depth/$S.unmapped.depth.txt" "$RUN/05_spades_output/contigs/$S.unmapped.fasta" > "$OUT/$S.high_copy.fa"
echo "$(grep -c '>' "$OUT/$S.high_copy.fa") contigs above ${MIN_DEPTH}x"

docker run --rm -v "$OUT":/o --user "$(id -u):$(id -g)" -e HOME=/tmp "$IMAGE" bash -c "
source /opt/conda/etc/profile.d/conda.sh; cd /o
conda activate tormes-1.3.0; prodigal -q -p meta -i $S.high_copy.fa -a $S.faa -o /dev/null; conda deactivate
conda activate humannenv
diamond blastp -q $S.faa -d /main/db/humann/uniref/uniref/uniref90_201901b_full.dmnd -o $S.diamond.tsv \
    -f 6 qseqid sseqid pident length qlen slen evalue bitscore -k 1 -e 1e-5 --threads $THREADS -b 2 -c 4 --tmpdir /tmp"
python3 "$HERE/name_hits.py" "$OUT/$S.diamond.tsv" > "$OUT/$S.named_hits.tsv"
echo "Wrote $OUT/$S.named_hits.tsv"
