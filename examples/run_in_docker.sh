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
#   --config  config template, path relative to the repository root. Its
#             step03.dir3 sets the host genome, which is required
#   --from/--to  first/last step (substring of the step script name)
#
# Environment variables:
#   IMAGE          Docker image (default metacheese)
#   INDEX_THREADS  threads for building the host Bowtie2 index (default 8)
#   STEP14_CHOICE  step 14 mode: 1 filter MAGs, 2 metadata for all, 3 both (default 1)
set -euo pipefail

IMAGE="${IMAGE:-metacheese}"
INDEX_THREADS="${INDEX_THREADS:-8}"
STEP14_CHOICE="${STEP14_CHOICE:-1}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

run_name="" input="" config=""
args=("$@")
for (( i=0; i<${#args[@]}; i++ )); do
    case "${args[$i]}" in
        --name)  run_name="${args[$((i+1))]}" ;;
        --input) input="${args[$((i+1))]}" ;;
        --config) config="${args[$((i+1))]}" ;;
    esac
done
[[ -n "$run_name" ]] || { echo "ERROR: --name is required" >&2; exit 1; }

docker image inspect "$IMAGE" >/dev/null 2>&1 \
    || { echo "ERROR: image '$IMAGE' not found. Run: docker build -t $IMAGE ." >&2; exit 1; }
if [[ -n "$input" ]]; then
    ls input/"$input"/*_1.fq.gz >/dev/null 2>&1 \
        || { echo "ERROR: no *_1.fq.gz files in input/$input" >&2; exit 1; }
fi
# The host comes from step03.dir3 of the config this run will use: the saved
# one when resuming, otherwise --config (same default as run_pipeline.sh).
if [[ -f "output/$run_name/config.yml" ]]; then
    config="output/$run_name/config.yml"
else
    config="${config:-config/Config.yml}"
fi
[[ -f "$config" ]] || { echo "ERROR: config not found: $config" >&2; exit 1; }
host_index="$(sed -n '/^step03:/,/^[^ ]/s/^[[:space:]]*dir3:[[:space:]]*"\{0,1\}\([^"#]*\).*/\1/p' "$config")"
host_index="${host_index%"${host_index##*[![:space:]]}"}"
# An absolute path is mounted at the same path in the container; paths under
# data/ are already mounted.
host_mount=()
case "$host_index" in
    *[[:space:]]*) echo "ERROR: Bowtie2 cannot use a host path with spaces: $host_index" >&2
       exit 1 ;;
    /*/*)        host_files="$host_index"; container_index="$host_index" ;;
    ../../data/*) host_files="${host_index#../../}"; container_index="/main/$host_files" ;;
    *) echo "ERROR: host removal needs a host genome: set step03.dir3 in $config to /path/to/host_folder/<index_name> (see README, Host genome)" >&2
       exit 1 ;;
esac
host_dir="${host_files%/*}"
if compgen -G "$host_files.1.bt2*" >/dev/null; then
    host_ro=":ro"   # nothing to build, so the host folder stays read-only
elif compgen -G "$host_dir/*.f*a" >/dev/null; then
    host_ro=""
else
    echo "ERROR: no host genome (.fa/.fna/.fasta) or Bowtie2 index $host_files in $host_dir (see README, Host genome)" >&2
    exit 1
fi
[[ "$host_index" != /* ]] || host_mount=(-v "$host_dir:$host_dir$host_ro")

mkdir -p "output/$run_name"
LOG="output/$run_name/run_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to $LOG"

# Same volume layout as docker-compose.yml. The container runs as root (some
# steps write into the conda environments); on exit the files it created are
# given back to the calling user.
docker run --rm \
    -e STEP14_CHOICE="$STEP14_CHOICE" \
    -e HOST_INDEX="$container_index" \
    -e INDEX_THREADS="$INDEX_THREADS" \
    -v "$REPO_ROOT/scripts:/main/scripts" \
    -v "$REPO_ROOT/config:/main/config" \
    -v "$REPO_ROOT/data:/main/data" \
    -v "$REPO_ROOT/docs:/main/docs" \
    -v "$REPO_ROOT/examples:/main/examples" \
    -v "$REPO_ROOT/input:/main/input" \
    -v "$REPO_ROOT/output:/main/output" \
    ${host_mount[@]+"${host_mount[@]}"} \
    -w /main \
    "$IMAGE" bash -c '
        set -euo pipefail
        trap "chown -R '"$(id -u):$(id -g)"' output/'"$run_name"' data/gene scripts/pipeline 2>/dev/null; find \"\${HOST_INDEX%/*}\" -maxdepth 1 -name \"\${HOST_INDEX##*/}.*.bt2*\" -user 0 -exec chown '"$(id -u):$(id -g)"' {} + 2>/dev/null || true" EXIT
        bash scripts/build_bowtie2_index.sh "$HOST_INDEX" "$INDEX_THREADS"
        bash scripts/run_pipeline.sh "$@"
    ' bash "$@" 2>&1 | tee "$LOG"
