#!/usr/bin/env bash
# Run the whole pipeline on the full PRJNA1198075 example (28 samples),
# unattended, inside the metacheese Docker image. Meant for a large server:
# check examples/PRJNA1198075/Config.yml (threads, SPAdes memory) first.
#
# Run from the repository root, on the host:
#   bash examples/PRJNA1198075/run_example.sh             # full run
#   bash examples/PRJNA1198075/run_example.sh --from 05   # resume from a step
#
# Prerequisites: see examples/PRJNA1198075/README.md
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

INDEX_THREADS="${INDEX_THREADS:-32}" exec bash examples/run_in_docker.sh \
    --input PRJNA1198075 \
    --name "${RUN_NAME:-example_PRJNA1198075}" \
    --config examples/PRJNA1198075/Config.yml \
    "$@"
