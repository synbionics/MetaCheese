#!/usr/bin/env bash
# Check a machine before running the full PRJNA1198075 example: Docker image,
# cores and memory against the config, free disk, input reads, host index.
# Changes nothing.
#
# Usage (repository root, on the host, after prepare_input.sh):
#   bash examples/PRJNA1198075/check_server.sh
# Environment: IMAGE (default metacheese), RUN_NAME (default example_PRJNA1198075)
#
# Exit code: 0 ready (maybe with warnings), 1 something must be fixed first.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
IMAGE="${IMAGE:-metacheese}"
RUN_NAME="${RUN_NAME:-example_PRJNA1198075}"
SAMPLES=28
TESTED=tests/regression/baseline/tables/versions.tsv   # tools and databases the regression test passed with
fails=0 warns=0
ok()   { printf 'OK    %s\n' "$*"; }
warn() { printf 'WARN  %s\n' "$*"; warns=$((warns + 1)); }
bad()  { printf 'FAIL  %s\n' "$*"; fails=$((fails + 1)); }

# A resumed run uses the config saved in its run folder, a new run the example's
CONFIG=examples/PRJNA1198075/Config.yml
if [[ -f "output/$RUN_NAME/config.yml" ]]; then
    CONFIG="output/$RUN_NAME/config.yml"
    ok "output/$RUN_NAME exists: the run resumes with its own $CONFIG"
    grep -q '^step17:' "$CONFIG" || bad "$CONFIG has no step17 section (made before the report became step 17): copy step16/step17 from examples/PRJNA1198075/Config.yml"
fi
cfg() {  # section key: value of "  key:" under "section:" in the config, comments removed
    awk -v s="$1:" -v k="$2:" '$0 ~ "^[^ #]" { in_s = ($1 == s) } in_s && $1 == k { $1 = ""; sub(/#.*/, ""); gsub(/[" ]/, ""); print; exit }' "$CONFIG"
}

echo "== Docker image ($IMAGE)"
if ! command -v docker >/dev/null; then
    bad "docker not found"
elif ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    bad "image $IMAGE not found: build it: docker build -t $IMAGE ."
else
    ok "image $IMAGE: $(docker image inspect -f '{{.Id}} created {{.Created}}' "$IMAGE" | cut -c8-19,72-)"
    # Tool and database versions against the tested ones (versions.sh prints
    # the same lines the regression test stores)
    found="$(docker run --rm -i "$IMAGE" bash -s < tests/regression/versions.sh 2>/dev/null)"
    if [[ -z "$found" ]]; then
        bad "could not read the tool versions of $IMAGE"
    else
        diffs="$(awk -F'\t' 'NR == FNR { t[$1] = $2; next } ($1 in t) && t[$1] != $2 { print "  " $1 ": " $2 " (tested: " t[$1] ")" }' \
            "$TESTED" - <<< "$found")"
        if [[ -z "$diffs" ]]; then
            ok "tool and database versions as tested ($(wc -l < "$TESTED") checked)"
        else
            warn "versions differ from the tested image:"; echo "$diffs"
        fi
    fi
fi

echo "== Cores and memory (config: $CONFIG)"
cores=$(nproc)
ram=$(awk '/^MemTotal:/ {printf "%d", $2 / 1048576}' /proc/meminfo)
max_threads=0
for pair in step00:par1 step00_01:par1 step03:par1 step04:par1 step04b_1:par1 step05:par1 step06:filter1 \
            step06b:par1 step07:par1 step08_09:par1 step10_11:par2 step12_13:par3 step15:par2; do
    t=$(cfg "${pair%%:*}" "${pair#*:}"); [[ "$t" =~ ^[0-9]+$ ]] && (( t > max_threads )) && max_threads=$t
done
if (( max_threads > cores )); then
    warn "config uses up to $max_threads threads, the machine has $cores cores: lower the par* threads"
else
    ok "up to $max_threads threads, $cores cores"
fi
spades=$(cfg step05 par2); block=$(cfg step04b_1 par2)
if [[ ! "$spades" =~ ^[0-9]+$ ]]; then
    bad "step05.par2 (metaSPAdes memory, GB) not found in $CONFIG"
elif (( spades >= ram )); then
    bad "metaSPAdes may use ${spades} GB (step05.par2), the machine has ${ram} GB: set it below the RAM"
else
    ok "metaSPAdes limit ${spades} GB of ${ram} GB RAM"
    (( spades < 100 )) && warn "the deepest samples (up to 85M read pairs) may need more than ${spades} GB in metaSPAdes"
fi
[[ "$block" =~ ^[0-9]+$ ]] && (( block * 6 > ram )) && warn "HUMAnN DIAMOND block size $block needs about $((block * 6)) GB, more than the ${ram} GB RAM"

echo "== Disk"
mkdir -p output
free_gb=$(df -BG --output=avail output | tail -1 | tr -dc 0-9)
if (( free_gb < 1000 )); then
    warn "${free_gb} GB free for output/: the full run needs about 1 TB"
else
    ok "${free_gb} GB free for output/"
fi

echo "== Input reads (input/PRJNA1198075)"
sheet=input/PRJNA1198075/samplesheet.tsv
if [[ ! -f "$sheet" ]]; then
    bad "$sheet not found: run examples/PRJNA1198075/prepare_input.sh <DATA_DIR>"
else
    n=0 missing=0
    while IFS=$'\t' read -r name _; do
        n=$((n + 1))
        for mate in 1 2; do
            [[ -s "input/PRJNA1198075/${name}_${mate}.fq.gz" ]] || { missing=$((missing + 1)); echo "      missing ${name}_${mate}.fq.gz"; }
        done
    done < <(tail -n +2 "$sheet")
    if (( n != SAMPLES || missing )); then
        bad "$n samples in $sheet (expected $SAMPLES), $missing FASTQ files missing or empty"
    else
        ok "$n samples, all FASTQ files present"
    fi
fi

echo "== Host genome and other data"
index=$(cfg step03 dir3)
case "$index" in
    ../../*) index_host="${index#../../}" ;;
    *) index_host="$index" ;;
esac
if compgen -G "$index_host.1.bt2*" >/dev/null; then
    ok "Bowtie2 host index $index_host"
elif compgen -G "${index_host%/*}/*.fna" >/dev/null || compgen -G "${index_host%/*}/*.fa" >/dev/null \
     || compgen -G "${index_host%/*}/*.fasta" >/dev/null; then
    ok "host FASTA in ${index_host%/*}/: the index is built at the start of the run (about 40 min)"
else
    bad "no host genome in ${index_host%/*}/: run bash examples/PRJNA1198075/download_host_genome.sh"
fi
[[ -s data/calculate_diversity.R ]] && ok "data/calculate_diversity.R" || bad "data/calculate_diversity.R missing (it is in git)"

echo "== Result: $fails to fix, $warns warnings"
(( fails == 0 ))
