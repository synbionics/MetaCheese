#!/usr/bin/env bash
# Run the whole pipeline on the small PRJNA1198075 example (6 samples x 5M
# read pairs), unattended, inside the metacheese Docker image.
#
# Run from the repository root, on the host:
#   bash examples/small_PRJNA1198075/run_example.sh                 # full run
#   bash examples/small_PRJNA1198075/run_example.sh --to 04.metaphlan
#   bash examples/small_PRJNA1198075/run_example.sh --from 05        # resume
#
# Prerequisites: see examples/small_PRJNA1198075/README.md
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

INDEX_THREADS="${INDEX_THREADS:-16}" exec bash examples/run_in_docker.sh \
    --input small_PRJNA1198075 \
    --name "${RUN_NAME:-example_small}" \
    --config examples/small_PRJNA1198075/Config.yml \
    "$@"
