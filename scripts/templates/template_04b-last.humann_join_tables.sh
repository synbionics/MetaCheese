#!/bin/bash
# Step 04b-last (diagram C.4, community-level aggregation): merge the
# per-sample HUMAnN outputs (C.3) into cohort tables, normalise them and
# regroup gene families into functions.
#
# Output (in the merged tables folder):
#   humann_joined_genefamilies.tsv        UniRef90 gene families, RPK, stratified by species
#   humann_joined_genefamilies_cpm.tsv    same, copies per million
#   humann_joined_pathabundance.tsv       MetaCyc pathways, RPK, stratified
#   humann_joined_pathabundance_cpm.tsv   same, copies per million
#   humann_joined_pathcoverage.tsv        MetaCyc pathway coverage (0-1)
#   merged_table_metacyc_rxn.tsv          gene families -> MetaCyc reactions (CPM, named)
#   merged_table_ko.tsv                   gene families -> KEGG Orthology (CPM, named) *
#   merged_table_ec.tsv                   gene families -> EC numbers (CPM, named) *
#   unstratified/*.tsv                    community totals of each table (no per-species rows)
#
# * needs the full HUMAnN utility mapping files (installed in the image since
#   the FastQC update). KEGG pathways are not produced: HUMAnN has no built-in
#   KO -> KEGG pathway mapping; use the MetaCyc pathway tables instead.
set -euo pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate humannenv

INPUT_DIR="@04b_last_dir1@"
OUTPUT_DIR="@04b_last_dir2@"

GENEFILES=("$INPUT_DIR"/*_genefamilies.tsv)
if [[ ! -e "${GENEFILES[0]}" ]]; then
    echo "Error: No *_genefamilies.tsv file found in $INPUT_DIR" >&2
    exit 2
fi

rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/unstratified"

# 1. Join per-sample tables into one table per output type
for kind in genefamilies pathabundance pathcoverage; do
    echo "Joining *_${kind}.tsv:"
    ls "$INPUT_DIR"/*_${kind}.tsv
    humann_join_tables \
        --input "$INPUT_DIR" \
        --file_name "_${kind}.tsv" \
        --output "$OUTPUT_DIR/humann_joined_${kind}.tsv"
done

# 2. Normalise abundances to copies per million (coverage is already 0-1)
for kind in genefamilies pathabundance; do
    humann_renorm_table \
        --input "$OUTPUT_DIR/humann_joined_${kind}.tsv" \
        --output "$OUTPUT_DIR/humann_joined_${kind}_cpm.tsv" \
        --units cpm --update-snames
done

# 3. Regroup gene families into functional categories and add readable names.
#    group = built-in regroup option, names = humann_rename_table option
# Read the help once: piping it into "grep -q" under pipefail can fail with
# SIGPIPE and wrongly report a grouping as missing.
REGROUP_HELP="$(humann_regroup_table --help 2>&1 || true)"
regroup() {
    local group="$1" names="$2" out="$3"
    if [[ "$REGROUP_HELP" != *"$group"* ]]; then
        echo "WARNING: grouping '$group' not available (HUMAnN utility_mapping files" \
             "not installed in this image), skipping $out" >&2
        return 0
    fi
    humann_regroup_table \
        --input "$OUTPUT_DIR/humann_joined_genefamilies_cpm.tsv" \
        --groups "$group" \
        --output "$OUTPUT_DIR/tmp_${out}"
    humann_rename_table \
        --input "$OUTPUT_DIR/tmp_${out}" \
        --names "$names" \
        --output "$OUTPUT_DIR/$out"
    rm -f "$OUTPUT_DIR/tmp_${out}"
}
regroup uniref90_rxn      metacyc-rxn    merged_table_metacyc_rxn.tsv
regroup uniref90_ko       kegg-orthology merged_table_ko.tsv
regroup uniref90_level4ec ec             merged_table_ec.tsv

# 4. Community-level (unstratified) copies, smaller and easier to plot
for table in "$OUTPUT_DIR"/*.tsv; do
    humann_split_stratified_table --input "$table" --output "$OUTPUT_DIR/unstratified" >/dev/null
done
rm -f "$OUTPUT_DIR"/unstratified/*_stratified.tsv
for f in "$OUTPUT_DIR"/unstratified/*_unstratified.tsv; do
    mv "$f" "${f%_unstratified.tsv}.tsv"
done

echo "Merged tables:"
ls -1 "$OUTPUT_DIR" "$OUTPUT_DIR/unstratified"
echo "HUMAnN merge and annotation completed."
conda deactivate
