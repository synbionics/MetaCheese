#!/bin/bash

# =============================
# 1. CONFIGURE WHAT TO DELETE/PRESERVE
# =============================

# List of subfolders (relative to the output folder) to process.
# All content will be deleted except what is specified in PRESERVE.
TO_DELETE=(
    "01_AdapterRemoval"
    #"03_bowtie2_output"
    "04_metaphlan_output"
    #"04b_humann_output"
    "05_spades_output"
    "07_Bowtie_Index"
    "08_mapping_coverage"
    #"10_metabat_depth"
    "11_metabat_MAG"
    "12_checkm"
    #"13_checkm2"
    #"14_MAGs_high_quality"
    #"15_tormes_MAGs"
)

# Map of exceptions (elements to preserve inside the respective folders).
# If you want to preserve the ENTIRE folder leave empty ("").
# Use commas for multiple elements/patterns. Supports wildcard/glob (e.g.: "*.2.gz").
declare -A PRESERVE

# Example: for 04_metaphlan_output preserve subfolder "diversity"
#PRESERVE["04_metaphlan_output"]="diversity"
# Example: for 03_bowtie2_output preserve all .2.gz files
#PRESERVE["03_bowtie2_output"]="*.2.gz"
# To preserve nothing, you can leave: PRESERVE["04b_humann_output"]=""
# You can add more entries (with ',' no spaces needed)

PRESERVE["04_metaphlan_output"]="diversity,merged_abundance_table.txt,*.txt"
PRESERVE["05_spades_output"]="contigs/filtered"

# =============================
# 2. SCRIPT LOGIC
# =============================

OUTPUT_BASE="output"

echo
ls -1d $OUTPUT_BASE/20* 2>/dev/null || { echo "No folder found in $OUTPUT_BASE!"; exit 1; }
echo
read -p "Output folder to clean (e.g.: 20250730_asdwq): " MAIN_DIR
MAIN_PATH="$OUTPUT_BASE/$MAIN_DIR"

if [[ ! -d "$MAIN_PATH" ]]; then
    echo "WARNING: Folder $MAIN_PATH does not exist."
    exit 2
fi

echo
echo "WARNING: The contents of the following folders will be DELETED:"
for FOLDER in "${TO_DELETE[@]}"; do
    PRES="${PRESERVE[$FOLDER]}"
    FOLDER_PATH="$MAIN_PATH/$FOLDER"
    if [[ -d "$FOLDER_PATH" ]]; then
        echo " - $FOLDER_PATH"
        if [[ -n "$PRES" ]]; then
            echo "   (will be preserved: $PRES)"
        fi
    else
        echo " - $FOLDER_PATH (NOT FOUND)"
    fi
done
echo

read -p "Do you confirm the selective deletion described above? [y/N]: " CONFIRM
CONFIRM=${CONFIRM,,}
if [[ "$CONFIRM" != "y" && "$CONFIRM" != "yes" ]]; then
    echo "Cancelled."
    exit 0
fi

# Function to delete everything except specified (pattern/glob supported)
clean_folder_except() {
    local DIR="$1"
    local PRESERVE_LIST="$2"
    [[ ! -d "$DIR" ]] && { echo "[NOT FOUND]: $DIR"; return; }

    # Build list of elements to preserve (also expands patterns)
    local EXCLUDES=()
    IFS=',' read -ra RAW_EXCLUDES <<< "$PRESERVE_LIST"
    for PATTERN in "${RAW_EXCLUDES[@]}"; do
        [[ -z "$PATTERN" ]] && continue
        for MATCH in "$DIR"/$PATTERN; do
            [[ -e "$MATCH" ]] && EXCLUDES+=("$MATCH")
        done
    done

    shopt -s dotglob nullglob
    for ITEM in "$DIR"/*; do
        SKIP=0
        for EX in "${EXCLUDES[@]}"; do
            [[ "$ITEM" == "$EX" ]] && SKIP=1 && break
        done
        if [[ $SKIP -eq 0 ]]; then
            rm -rf "$ITEM"
            echo "Deleted: $ITEM"
        else
            echo "Preserved: $ITEM"
        fi
    done
    shopt -u dotglob nullglob
}

# Actual cleaning execution
for FOLDER in "${TO_DELETE[@]}"; do
    DIR="$MAIN_PATH/$FOLDER"
    PRES="${PRESERVE[$FOLDER]}"
    if [[ -d "$DIR" ]]; then
        if [[ -n "$PRES" ]]; then
            clean_folder_except "$DIR" "$PRES"
        else
            rm -rf "$DIR"
            echo "Deleted EVERYTHING: $DIR"
        fi
    else
        echo "[NOT FOUND]: $DIR"
    fi
done

echo
echo "Cleaning completed!"
