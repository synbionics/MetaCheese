#!/usr/bin/env bash
# Build the small MetaCheese example: 3 PDO and 3 non-PDO cheese samples from
# PRJNA1198075, first N read pairs of each, written into
# input/small_PRJNA1198075/ with the file naming the pipeline expects
# (<sample>_1.fq.gz / <sample>_2.fq.gz) plus a samplesheet.tsv.
#
# By default the reads are streamed from ENA and only the first N pairs of the
# 12 FASTQ files are downloaded (about 3.5 GB for 5M pairs), so the full
# dataset (125 GB) is not needed. If you already downloaded it with
# examples/PRJNA1198075/download_PRJNA1198075.sh, pass that folder to read the
# local files instead.
#
# Keeping the first N reads (instead of random subsampling) is fast and fully
# deterministic: both modes give the same reads. Reads at the start of a FASTQ
# come from the first flow-cell tiles, which can be slightly lower quality,
# but this does not bias the taxonomic composition in any meaningful way.
#
# Usage (from the repository root, on the host - no container needed):
#   bash examples/small_PRJNA1198075/prepare_input.sh               # stream from ENA
#   bash examples/small_PRJNA1198075/prepare_input.sh <FULL_DATA_DIR>   # local files
#
# Environment variables:
#   READ_PAIRS  read pairs to keep per sample (default 5000000)
#   OUT_DIR     output folder (default input/small_PRJNA1198075)
set -euo pipefail

SRC="${1:-}"
PAIRS="${READ_PAIRS:-5000000}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${OUT_DIR:-$REPO_ROOT/input/small_PRJNA1198075}"

# 3 PDO + 3 non-PDO samples, as chosen in the METACHEESE manuscript (dataset
# example): D_GP_04 + ND_GP_4 and D_GP_02 + ND_GP_10 are well separated in the
# beta-diversity plot, while D_GP_13 + ND_GP_15 are similar to each other.
# sample_name  run_accession  biosample  title  group (from ENA + NCBI BioSample)
SAMPLES="\
D_GP_02  SRR31729973 SAMN45828310 cheese02  PDO
D_GP_04  SRR31729955 SAMN45828335 cheese04  PDO
D_GP_13  SRR31729948 SAMN45828317 cheese13  PDO
ND_GP_4  SRR31729957 SAMN45828333 Ncheese4  NONPDO
ND_GP_10 SRR31729968 SAMN45828323 Ncheese10 NONPDO
ND_GP_15 SRR31729961 SAMN45828329 Ncheese15 NONPDO"

if command -v pigz >/dev/null; then GZ="pigz -p 4"; ZCAT="pigz -dc"; else GZ="gzip"; ZCAT="zcat"; fi
LINES=$((PAIRS * 4))

if [[ -n "$SRC" ]]; then
    [[ -d "$SRC" ]] || { echo "ERROR: $SRC is not a folder" >&2; exit 1; }
    echo "Reading local files from $SRC"
else
    echo "Streaming the first $PAIRS read pairs of each sample from ENA"
fi
mkdir -p "$OUT"
echo "Writing $PAIRS read pairs per sample to $OUT"

# Write the first $LINES lines of a gzipped FASTQ (local path or URL) to $2.
# head closes the pipe early on purpose, so the SIGPIPE/curl write error of
# the producer is ignored; the line count check below catches real failures
# (missing file, dropped connection).
take_head() {
    local from="$1" to="$2"
    if [[ "$from" == https://* ]]; then
        { curl -sSfL --retry 5 "$from" 2>/dev/null || true; } | { $ZCAT 2>/dev/null || true; } \
            | head -n "$LINES" | $GZ > "$to"
    else
        { $ZCAT "$from" || true; } | head -n "$LINES" | $GZ > "$to"
    fi
}

printf 'sample_name\trun_accession\tbiosample\ttitle\tgroup\tfastq_1\tfastq_2\tread_pairs\n' > "$OUT/samplesheet.tsv.tmp"

while read -r s run biosample title group; do
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$run" "$biosample" "$title" "$group" \
        "${s}_1.fq.gz" "${s}_2.fq.gz" "$PAIRS" >> "$OUT/samplesheet.tsv.tmp"

    urls=()
    if [[ -z "$SRC" ]]; then
        # FASTQ URLs of the run, as listed by ENA (mate 1 and mate 2)
        ftp="$(curl -sSf --retry 5 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp&format=tsv" \
            | awk -F'\t' 'NR==2 {print $2}')"
        IFS=';' read -r -a urls <<< "$ftp"
        (( ${#urls[@]} == 2 )) || { echo "ERROR: ENA did not return two FASTQ files for $run" >&2; exit 1; }
    fi

    for mate in 1 2; do
        dst_file="$OUT/${s}_${mate}.fq.gz"
        if [[ -s "$dst_file" ]] && gzip -t "$dst_file" 2>/dev/null \
            && (( $($ZCAT "$dst_file" | wc -l) == LINES )); then
            echo "  $dst_file exists, skipping"
            continue
        fi
        if [[ -n "$SRC" ]]; then
            from="$SRC/${run}_${mate}.fastq.gz"
            [[ -s "$from" ]] || { echo "ERROR: missing $from (re-run the download script)" >&2; exit 1; }
        else
            from="https://${urls[mate - 1]}"
        fi
        echo "  $s mate $mate <- $(basename "$from")"

        ok=0
        for attempt in 1 2 3; do
            take_head "$from" "$dst_file.tmp"
            got=$($ZCAT "$dst_file.tmp" | wc -l)
            if (( got == LINES )); then ok=1; break; fi
            echo "  got $((got / 4)) reads instead of $PAIRS (attempt $attempt)" >&2
            [[ -n "$SRC" ]] && break
            sleep 30
        done
        (( ok )) || { rm -f "$dst_file.tmp"; echo "ERROR: could not get $PAIRS read pairs for $s mate $mate" >&2; exit 1; }
        mv "$dst_file.tmp" "$dst_file"
    done
done <<< "$SAMPLES"

mv "$OUT/samplesheet.tsv.tmp" "$OUT/samplesheet.tsv"
echo "Done. Sample sheet: $OUT/samplesheet.tsv"
