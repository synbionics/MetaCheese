#!/usr/bin/env bash
# Screen every sample of a dataset cheaply, before choosing samples or spending
# compute: a high-copy k-mer proxy, and MetaPhlAn on the first read pairs.
# Took about 2 h for 28 samples with 8 threads.
#
# Usage: bash experiments/branch_gain/screen_samples.sh <DATA_DIR> <outdir>
#   <DATA_DIR> holds samplesheet.tsv and the FASTQ files (examples/PRJNA1198075)
# Environment: PAIRS (default 1000000), KMER_READS (default 500000), THREADS (default 8)
set -euo pipefail
D="$(realpath "$1")"; OUT="$2"; mkdir -p "$OUT/mpa"; OUT="$(realpath "$OUT")"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAIRS="${PAIRS:-1000000}"; KMER_READS="${KMER_READS:-500000}"; THREADS="${THREADS:-8}"

export PYTHONHASHSEED=0
tail -n +2 "$D/samplesheet.tsv" | awk -F'\t' '{print $1, $6}' \
    | xargs -P "$THREADS" -L 1 bash -c 'python3 -I "'"$HERE"'/kmer_highcopy.py" "$0" "'"$D"'/$1" '"$KMER_READS" > "$OUT/kmer.tsv"

docker run --rm -v "$D":/d:ro -v "$OUT":/o --user "$(id -u):$(id -g)" -e HOME=/tmp metacheese bash -c '
source /opt/conda/etc/profile.d/conda.sh; conda activate bioenv
tail -n +2 /d/samplesheet.tsv | while IFS=$'"'"'\t'"'"' read -r n run bs t g f1 f2 p; do
    [ -s /o/mpa/$n.txt ] && continue
    zcat /d/$f1 2>/dev/null | head -n '"$((PAIRS * 4))"' > /tmp/$n.1.fq || true
    zcat /d/$f2 2>/dev/null | head -n '"$((PAIRS * 4))"' > /tmp/$n.2.fq || true
    metaphlan /tmp/$n.1.fq,/tmp/$n.2.fq --input_type fastq --bowtie2db /main/db/metaphlan \
        --index mpa_vJun23_CHOCOPhlAnSGB_202307 --nproc '"$THREADS"' --bowtie2out /tmp/$n.bt2.bz2 -o /o/mpa/$n.txt > /dev/null 2>> /o/mpa.err
    rm -f /tmp/$n.*
done
merge_metaphlan_tables.py /o/mpa/*.txt > /o/merged.txt'
python3 "$HERE/screen_table.py" "$OUT/merged.txt" "$OUT/kmer.tsv" "$D/samplesheet.tsv" | tee "$OUT/screen_table.tsv"
