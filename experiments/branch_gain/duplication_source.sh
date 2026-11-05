#!/usr/bin/env bash
# Are a sample's duplicated reads real high-copy DNA or a library (PCR) artefact?
#
# Counts read pairs with identical contig, position, strand and insert size,
# per class of contig depth. Library duplication affects contigs of every
# depth; real high-copy DNA gives duplicates only on the deep contigs, because
# those contigs run out of distinct start positions.
#
# Usage (repository root, on the host):
#   bash experiments/branch_gain/duplication_source.sh <run_folder> <sample> [<sample> ...]
set -euo pipefail
RUN="$(realpath "$1")"; shift
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
docker run --rm -i -v "$RUN":/r:ro -v "$HERE":/x:ro metacheese bash -c '
source /opt/conda/etc/profile.d/conda.sh; conda activate bioenv
for s in "$@"; do
    samtools view -f 66 -F 2316 /r/08_mapping_coverage/$s.unmapped.sorted.bam \
        | python3 /x/duplication_count.py /r/10_metabat_depth/$s.unmapped.depth.txt "$s"
done' bash "$@" 2>&1 | grep -v ttyname
