#!/usr/bin/env bash
# Step 06b (diagram A.8): quality of the filtered assemblies with QUAST.
# All samples go into one QUAST run, so report.html / report.tsv compare them
# side by side. No reference genome is given (metagenome), so QUAST reports
# the reference-free statistics (contigs, total length, N50/L50, largest
# contig, GC) and downloads nothing at run time.
set -euo pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate quastenv

contig_dir="@06b_dir1@"
output_dir="@06b_dir2@"
mkdir -p "$output_dir"

shopt -s nullglob
contigs=("$contig_dir"/*.fasta_sort.fasta)
if (( ${#contigs[@]} == 0 )); then
    echo "ERROR: no filtered contigs (*.fasta_sort.fasta) in $contig_dir" >&2
    exit 1
fi

# Sample names as labels (D_GP_02.unmapped.fasta_sort.fasta -> D_GP_02)
labels=$(for f in "${contigs[@]}"; do
    n=$(basename "$f" .fasta_sort.fasta); echo "${n%.unmapped}"
done | paste -sd, -)

echo "[$(date)] QUAST on ${#contigs[@]} assemblies -> $output_dir"
quast.py "${contigs[@]}" \
    --labels "$labels" \
    --output-dir "$output_dir" \
    --threads @06b_par1@ \
    --min-contig @06b_par2@

echo "[$(date)] QUAST completed: $output_dir/report.html"
conda deactivate
