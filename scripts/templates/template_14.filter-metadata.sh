#!/bin/bash
set -euo pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate bioenv

# Common parameters
TSV_FILE="@14_var1@"
SOURCE_DIR="@14_var2@"
TARGET_DIR="@14_var3@"
ALL_MAGS_DIR="@14_var4@"
COMPLETENESS="@14_var5@"
CONTAMINATION="@14_var6@"
METADATA_FILTERED="@14_var7@"
METADATA_ALL="@14_var8@"

# Input checks
[[ -f "$TSV_FILE" ]] || { echo "Error: TSV file $TSV_FILE not found"; exit 1; }
[[ -d "$SOURCE_DIR" ]] || { echo "Error: SOURCE_DIR $SOURCE_DIR not found"; exit 1; }

# The choice can be given as first argument or via STEP14_CHOICE (used by the
# non-interactive runner); otherwise it is asked interactively.
scelta="${1:-${STEP14_CHOICE:-}}"
if [[ -z "$scelta" ]]; then
    echo "Choose the process:"
    echo "1) Filter MAGs by quality and generate metadata only for filtered ones"
    echo "2) Generate metadata for all MAGs in a directory (no filter)"
    echo "3) Run both"
    read -p "Choice (1/2/3): " scelta
fi

if [[ "$scelta" != "1" && "$scelta" != "2" && "$scelta" != "3" ]]; then
    echo "Error: invalid choice" >&2
    exit 1
fi

if [[ "$scelta" == "1" || "$scelta" == "3" ]]; then
    echo "Running filter and generating filtered metadata..."
    mkdir -p "$TARGET_DIR"
    FILTERED_TSV="filtered_data.tsv"

    awk -F'\t' "NR==1 || (\$2 > $COMPLETENESS && \$3 < $CONTAMINATION)" "$TSV_FILE" > "$FILTERED_TSV"

    if [[ $(wc -l < "$FILTERED_TSV") -le 1 ]]; then
        echo "Error: no MAGs passed the filters (completeness > $COMPLETENESS, contamination < $CONTAMINATION)." >&2
        exit 1
    fi

    awk -F'\t' 'NR > 1 {print $1".fa"}' "$FILTERED_TSV" | while IFS= read -r filename; do
        if [[ -e "$SOURCE_DIR/$filename" ]]; then
            mv "$SOURCE_DIR/$filename" "$TARGET_DIR/"
            echo "File $filename moved to $TARGET_DIR"
        else
            echo "WARNING: file $filename not found in $SOURCE_DIR" >&2
        fi
    done

    ls "$TARGET_DIR"/*.fa > /dev/null || { echo "Error: no .fa file moved to $TARGET_DIR" >&2; exit 1; }

    ls "$TARGET_DIR"/*.fa | sed "s|.*/||; s/.fa//" > uno.tmp
    while read -r filename; do
        echo "GENOME" >> dos.tmp
        realpath "$TARGET_DIR/${filename}.fa" >> tres.tmp
        echo "This is the genome for ${filename}" >> cuatro.tmp
    done < uno.tmp
    paste uno.tmp dos.tmp tres.tmp cuatro.tmp | sed "1iSamples\tRead1\tRead2\tDescription" > "$TARGET_DIR/$METADATA_FILTERED"
    rm uno.tmp dos.tmp tres.tmp cuatro.tmp
    echo "Filtering completed: metadata written to $TARGET_DIR/$METADATA_FILTERED"
fi

if [[ "$scelta" == "2" || "$scelta" == "3" ]]; then
    echo "Generating metadata for all MAGs..."
    ls "$ALL_MAGS_DIR"/*.fa > /dev/null || { echo "Error: no .fa file found in $ALL_MAGS_DIR" >&2; exit 1; }

    ls "$ALL_MAGS_DIR"/*.fa | sed "s|.*/||; s/.fa//" > uno.tmp
    while read -r filename; do
        echo "GENOME" >> dos.tmp
        realpath "$ALL_MAGS_DIR/${filename}.fa" >> tres.tmp
        echo "This is the genome for ${filename}" >> cuatro.tmp
    done < uno.tmp
    paste uno.tmp dos.tmp tres.tmp cuatro.tmp | sed "1iSamples\tRead1\tRead2\tDescription" > "$ALL_MAGS_DIR/$METADATA_ALL"
    rm uno.tmp dos.tmp tres.tmp cuatro.tmp
    echo "Metadata for all MAGs written to $ALL_MAGS_DIR/$METADATA_ALL"
fi

echo "Process completed successfully."
conda deactivate
