#!/usr/bin/env bash
# Build the small MetaCheese example from the full PRJNA1198075 dataset.
#
# Takes 3 PDO and 3 non-PDO cheese samples and keeps the first N read pairs of
# each, writing them into input/<example_name>/ with the file naming the
# pipeline expects (<sample>_1.fq.gz / <sample>_2.fq.gz).
#
# Keeping the first N reads (instead of random subsampling) is fast and fully
# deterministic. Reads at the start of a FASTQ come from the first flow-cell
# tiles, which can be slightly lower quality, but this does not bias the
# taxonomic composition in any meaningful way.
#
# Usage (from the repository root, on the host - no container needed):
#   bash examples/small_PRJNA1198075/prepare_input.sh <FULL_DATA_DIR> [READ_PAIRS]
#
#   FULL_DATA_DIR  folder with the downloaded PRJNA1198075 FASTQ files and
#                  samplesheet.tsv (see examples/PRJNA1198075/download_PRJNA1198075.sh)
#   READ_PAIRS     read pairs to keep per sample (default 5000000)
set -euo pipefail

SRC="${1:?Usage: $0 <FULL_DATA_DIR> [READ_PAIRS]}"
PAIRS="${2:-5000000}"
EXAMPLE_NAME="small_PRJNA1198075"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO_ROOT/input/$EXAMPLE_NAME"

# 3 PDO + 3 non-PDO samples, as chosen in the METACHEESE manuscript (dataset
# example): D_GP_04 + ND_GP_4 and D_GP_02 + ND_GP_10 are well separated in the
# beta-diversity plot, while D_GP_13 + ND_GP_15 are similar to each other.
SAMPLES=(D_GP_02 D_GP_04 D_GP_13 ND_GP_4 ND_GP_10 ND_GP_15)

SHEET="$SRC/samplesheet.tsv"
[[ -f "$SHEET" ]] || { echo "ERROR: $SHEET not found" >&2; exit 1; }

if command -v pigz >/dev/null; then GZ="pigz -p 4"; ZCAT="pigz -dc"; else GZ="gzip"; ZCAT="zcat"; fi
LINES=$((PAIRS * 4))

mkdir -p "$OUT"
echo "Writing $PAIRS read pairs per sample to $OUT"

# Sample sheet for this example, same columns as the full one
head -n1 "$SHEET" > "$OUT/samplesheet.tsv"

for s in "${SAMPLES[@]}"; do
    row="$(awk -F'\t' -v s="$s" '$1==s' "$SHEET")"
    [[ -n "$row" ]] || { echo "ERROR: sample $s not in $SHEET" >&2; exit 1; }
    printf '%s\n' "$row" | awk -F'\t' -v OFS='\t' -v n="$PAIRS" \
        '{$6=$1"_1.fq.gz"; $7=$1"_2.fq.gz"; $8=n; print}' >> "$OUT/samplesheet.tsv"

    for mate in 1 2; do
        src_file="$SRC/$(printf '%s\n' "$row" | cut -f$((5 + mate)))"
        dst_file="$OUT/${s}_${mate}.fq.gz"
        [[ -s "$src_file" ]] || { echo "ERROR: missing source file $src_file (re-run the download script)" >&2; exit 1; }
        if [[ -s "$dst_file" ]] && gzip -t "$dst_file" 2>/dev/null; then
            echo "  $dst_file exists, skipping"
            continue
        fi
        echo "  $s mate $mate <- $(basename "$src_file")"
        # head closes the pipe early on purpose; ignore the resulting SIGPIPE
        { $ZCAT "$src_file" || true; } | head -n "$LINES" | $GZ > "$dst_file.tmp"
        got=$($ZCAT "$dst_file.tmp" | wc -l)
        (( got == LINES )) || { echo "ERROR: $dst_file has $((got / 4)) reads, expected $PAIRS" >&2; exit 1; }
        mv "$dst_file.tmp" "$dst_file"
    done
done

echo "Done. Sample sheet: $OUT/samplesheet.tsv"
