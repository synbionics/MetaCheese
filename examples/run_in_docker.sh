#!/usr/bin/env bash
# Run the MetaCheese pipeline unattended inside the metacheese Docker image.
# Used by the examples, but works for any input folder.
#
# Run from the host:
#   bash examples/run_in_docker.sh --input <input_folder> --name <run_name> \
#        --config <Config.yml> [--from <step>] [--to <step>]
#
#   --input   folder under input/ with <sample>_1.fq.gz / <sample>_2.fq.gz
#   --name    run folder under output/ (re-use the same name to resume)
#   --config  config template, path relative to the repository root
#   --from/--to  first/last step (substring of the step script name)
#
# Environment variables:
#   IMAGE          Docker image (default metacheese)
#   HOST_SPECIES   folder under data/gene/ with the host genome (default Bos_taurus)
#   INDEX_THREADS  threads for building the host Bowtie2 index (default 8)
#   STEP14_CHOICE  step 14 mode: 1 filter MAGs, 2 metadata for all, 3 both (default 1)
set -euo pipefail

IMAGE="${IMAGE:-metacheese}"
HOST_SPECIES="${HOST_SPECIES:-Bos_taurus}"
INDEX_THREADS="${INDEX_THREADS:-8}"
STEP14_CHOICE="${STEP14_CHOICE:-1}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

run_name="" input=""
args=("$@")
for (( i=0; i<${#args[@]}; i++ )); do
    case "${args[$i]}" in
        --name)  run_name="${args[$((i+1))]}" ;;
        --input) input="${args[$((i+1))]}" ;;
    esac
done
[[ -n "$run_name" ]] || { echo "ERROR: --name is required" >&2; exit 1; }

docker image inspect "$IMAGE" >/dev/null 2>&1 \
    || { echo "ERROR: image '$IMAGE' not found. Run: docker build -t $IMAGE ." >&2; exit 1; }
if [[ -n "$input" ]]; then
    ls input/"$input"/*_1.fq.gz >/dev/null 2>&1 \
        || { echo "ERROR: no *_1.fq.gz files in input/$input" >&2; exit 1; }
fi
host_dir="data/gene/$HOST_SPECIES"
compgen -G "$host_dir/*.f*a" >/dev/null || [[ -f "$host_dir/$HOST_SPECIES.1.bt2" ]] \
    || { echo "ERROR: no host genome (.fa/.fna/.fasta) in $host_dir (see README, Host genome)" >&2; exit 1; }

mkdir -p "output/$run_name"
LOG="output/$run_name/run_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to $LOG"

# Same volume layout as docker-compose.yml. The container runs as root (some
# steps write into the conda environments); on exit the files it created are
# given back to the calling user.
docker run --rm \
    -e STEP14_CHOICE="$STEP14_CHOICE" \
    -v "$REPO_ROOT/scripts:/main/scripts" \
    -v "$REPO_ROOT/config:/main/config" \
    -v "$REPO_ROOT/data:/main/data" \
    -v "$REPO_ROOT/examples:/main/examples" \
    -v "$REPO_ROOT/input:/main/input" \
    -v "$REPO_ROOT/output:/main/output" \
    -w /main \
    "$IMAGE" bash -c '
        set -euo pipefail
        trap "chown -R '"$(id -u):$(id -g)"' output/'"$run_name"' data/gene scripts/pipeline 2>/dev/null || true" EXIT
        bash scripts/build_bowtie2_index.sh '"$HOST_SPECIES $INDEX_THREADS"'
        bash scripts/run_pipeline.sh "$@"
    ' bash "$@" 2>&1 | tee "$LOG"
