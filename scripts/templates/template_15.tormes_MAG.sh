#!/bin/bash
set -eo pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate tormes-1.3.0

WORKDIR="@15_var1@"
INPUT="@15_var2@"
OUTPUT="@15_par1@"
THREADS="@15_par2@"

cd "$WORKDIR"

[[ -d "$WORKDIR" ]] || { echo "Error: working directory $WORKDIR not found" >&2; exit 1; }
[[ -f "$INPUT" ]] || { echo "Error: metadata file $INPUT not found" >&2; exit 1; }

echo "Running TORMES on $INPUT in directory $WORKDIR, output to $OUTPUT using $THREADS threads"

tormes --metadata "$INPUT" --output "$OUTPUT" --threads "$THREADS"

echo "TORMES successfully completed in $WORKDIR/$OUTPUT"

# TORMES 1.3.0 cannot build its summary report (tormes_report.html) when the MAGs
# have MLST schemes with different numbers of loci: R reads mlst.tab with as many
# columns as its first lines, so longer rows wrap, and two wrapped rows with the
# same first locus stop the report. Pad the report's copy of mlst.tab (TORMES'
# own mlst/mlst.tab is not changed) to one column count and render it again.
REPORT_FILES="$OUTPUT/report_files"
if [[ ! -s "$OUTPUT/tormes_report.html" && -s "$REPORT_FILES/mlst.tab" && -f "$REPORT_FILES/render_report.sh" ]]; then
    echo "TORMES summary report missing: padding $REPORT_FILES/mlst.tab and rendering it again"
    awk -F'\t' -v OFS='\t' 'NR == FNR { if (NF > n) n = NF; next } { $n = $n; print }' \
        "$REPORT_FILES/mlst.tab" "$REPORT_FILES/mlst.tab" > "$REPORT_FILES/mlst.tab.tmp"
    mv "$REPORT_FILES/mlst.tab.tmp" "$REPORT_FILES/mlst.tab"
    if bash "$REPORT_FILES/render_report.sh" > "$REPORT_FILES/render_report.log" 2>&1 \
        && [[ -s "$REPORT_FILES/tormes_report.html" ]]; then
        cp -f "$REPORT_FILES/tormes_report.html" "$OUTPUT/tormes_report.html"
        echo "TORMES summary report written to $OUTPUT/tormes_report.html"
    else
        echo "WARNING: the TORMES summary report could still not be built (see $REPORT_FILES/render_report.log);" \
             "the per-MAG TORMES results are complete" >&2
    fi
fi

conda deactivate
