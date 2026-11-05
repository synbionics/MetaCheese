#!/usr/bin/env bash
# Settle species names independently of both branches: compare every bin (kept
# or rejected) with an NCBI reference genome of each species that MetaPhlAn
# reports at >= MIN_ABUND % in the run. Mash is in the image and gives an estimate;
# confirm the values that matter, especially near 95%, with fastANI.
#
# Usage (repository root, on the host; needs internet):
#   bash experiments/branch_gain/species_check.sh <run_folder> <outdir>
# Environment: MIN_ABUND (default 1), IMAGE (default metacheese)
set -euo pipefail
RUN="$(realpath "$1")"; OUT="$2"; mkdir -p "$OUT/refs"; OUT="$(realpath "$OUT")"
MIN_ABUND="${MIN_ABUND:-1}"; IMAGE="${IMAGE:-metacheese}"
API=https://api.ncbi.nlm.nih.gov/datasets/v2/genome

awk -F'\t' -v m="$MIN_ABUND" '!/^#/ && $1 ~ /s__/ && $1 !~ /t__/ {
    for (i = 2; i <= NF; i++) if ($i + 0 >= m) { n = $1; sub(/.*s__/, "", n); print n; break } }' \
    "$RUN/04_metaphlan_output/merged_abundance_table.txt" | sort -u > "$OUT/species.txt"

while read -r sp; do
    fa="$OUT/refs/$sp.fna"; [[ -s "$fa" ]] && continue
    acc=$(curl -sS -m 30 "$API/taxon/${sp//_/%20}/dataset_report?filters.reference_only=true&page_size=1" \
        | python3 -c 'import json,sys; print(json.load(sys.stdin)["reports"][0]["accession"])' 2>/dev/null || true)
    if [[ -z "$acc" ]]; then echo "no reference genome: $sp" | tee -a "$OUT/missing_refs.txt"; continue; fi
    tmp="$(mktemp -d "$OUT/refs/tmp.XXXX")"
    curl -sS -m 300 -o "$tmp/g.zip" "$API/accession/$acc/download?include_annotation_type=GENOME_FASTA"
    unzip -q -j "$tmp/g.zip" '*_genomic.fna' -d "$tmp"
    cat "$tmp"/*_genomic.fna > "$fa"; rm -rf "$tmp"
    printf '%s\t%s\n' "$sp" "$acc" >> "$OUT/reference_accessions.tsv"
done < "$OUT/species.txt"

docker run --rm -v "$RUN":/r:ro -v "$OUT":/o --user "$(id -u):$(id -g)" -e HOME=/tmp "$IMAGE" bash -c '
source /opt/conda/etc/profile.d/conda.sh; conda activate tormes-1.3.0; cd /o
mash sketch -s 10000 -k 21 -p 4 -o refs refs/*.fna 2>/dev/null
for d in 14_MAGs_high_quality 11_metabat_MAG; do
    for b in /r/$d/*.fa; do
        [ -e "$b" ] || continue
        n=$(basename "$b" .fa); [ $d = 14_MAGs_high_quality ] && k=kept || k=rejected
        mash dist -p 2 -s 10000 refs.msh "$b" 2>/dev/null \
            | awk -v n="$n" -v k="$k" "{f=\$1; sub(/.*\//,\"\",f); sub(/.fna\$/,\"\",f); printf \"%s\t%s\t%s\t%.2f\t%s\n\",n,k,f,100*(1-\$3),\$5}" \
            | sort -t"	" -k4,4nr | head -2
    done
done' 2>&1 | grep -v ttyname > "$OUT/.mash.tsv"
{ printf 'bin\tstatus\treference\tmash_ani_pct\tshared_hashes\n'; cat "$OUT/.mash.tsv"; } > "$OUT/bins_vs_references.tsv"
rm -f "$OUT/.mash.tsv"
echo "Wrote $OUT/bins_vs_references.tsv (two best references per bin)"
