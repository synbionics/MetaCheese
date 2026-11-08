#!/usr/bin/env bash
# Run the whole workflow on the regression samples, summarise the results and
# compare them with the stored baseline, if there is one.
#
# Usage (repository root, on the host):
#   bash tests/regression/run.sh <run_name> [--from <step>] [--to <step>]
#
# Results: output/<run_name>/test_summary/ and output/<run_name>/test_comparison.md
#
# Environment variables: IMAGE (default metacheese), INDEX_THREADS (default 16)
#
# Exit code: 0 pass, 1 failed, 2 differences to check by hand.
set -uo pipefail

run_name="${1:?Usage: $0 <run_name> [--from <step>] [--to <step>]}"; shift
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
IMAGE="${IMAGE:-metacheese}"
TESTS=tests/regression

[[ -f input/test_regression/samplesheet.tsv ]] \
    || { echo "ERROR: no test input, run $TESTS/prepare_input.sh first" >&2; exit 1; }

IMAGE="$IMAGE" INDEX_THREADS="${INDEX_THREADS:-16}" bash examples/run_in_docker.sh \
    --input test_regression --name "$run_name" \
    --config examples/small_PRJNA1198075/Config.yml "$@"
pipeline_exit=$?
[[ -d "output/$run_name" ]] || exit 1

docker run --rm -i "$IMAGE" bash -s < "$TESTS/versions.sh" > "output/$run_name/versions.tsv"
python3 "$TESTS/summarize.py" "output/$run_name" --pipeline-exit "$pipeline_exit" || exit 1

if [[ -f "$TESTS/baseline/summary.json" ]]; then
    python3 "$TESTS/compare.py" "$TESTS/baseline/summary.json" \
        "output/$run_name/test_summary/summary.json" --out "output/$run_name/test_comparison.md"
    exit $?
fi
echo "No baseline stored yet. To make this run the baseline:"
echo "  rm -rf $TESTS/baseline && cp -r output/$run_name/test_summary $TESTS/baseline"
exit "$pipeline_exit"
