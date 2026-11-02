#!/usr/bin/env bash
# Prepare the full PRJNA1198075 example (28 cheese samples) as pipeline input.
#
# Files are placed in input/PRJNA1198075/ named <sample>_1.fq.gz / _2.fq.gz,
# as the pipeline expects. Hard links are used when the data is on the same
# filesystem as the repository (no extra disk space), otherwise files are copied.
#
# Usage (from the repository root, on the host):
#   bash examples/PRJNA1198075/download_PRJNA1198075.sh PRJNA1198075 <DATA_DIR>
#   bash examples/PRJNA1198075/prepare_input.sh <DATA_DIR>
set -euo pipefail

SRC="${1:?Usage: $0 <DATA_DIR>   (folder created by download_PRJNA1198075.sh)}"
SRC="$(cd "$SRC" && pwd)"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$REPO_ROOT/input/PRJNA1198075"
SHEET="$SRC/samplesheet.tsv"
[[ -f "$SHEET" ]] || { echo "ERROR: $SHEET not found" >&2; exit 1; }

mkdir -p "$OUT"
head -n1 "$SHEET" > "$OUT/samplesheet.tsv"
missing=0
while IFS=$'\t' read -r name run biosample title group fq1 fq2 pairs; do
    for mate in 1 2; do
        src="$SRC/$([[ $mate == 1 ]] && echo "$fq1" || echo "$fq2")"
        dst="$OUT/${name}_${mate}.fq.gz"
        if [[ ! -s "$src" ]]; then
            echo "MISSING $src" >&2; missing=1; continue
        fi
        [[ -e "$dst" ]] && continue
        ln "$src" "$dst" 2>/dev/null || cp "$src" "$dst"
    done
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$run" "$biosample" "$title" "$group" \
        "${name}_1.fq.gz" "${name}_2.fq.gz" "$pairs" >> "$OUT/samplesheet.tsv"
done < <(tail -n +2 "$SHEET")

(( missing == 0 )) || { echo "ERROR: some files are missing; re-run the download script (it resumes)" >&2; exit 1; }
echo "Done: $(( $(wc -l < "$OUT/samplesheet.tsv") - 1 )) samples in $OUT"
