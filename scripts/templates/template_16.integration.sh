#!/bin/bash
# Step R.1 (joint step of both branches): analyses that combine the results of the
# community-level branch (C.1-C.4) and of the assembly-specific branch (A.1-A.7).
# Writes tables to <run>/16_integration/ (read_fate.tsv, mags.tsv, notes.txt),
# read by the report (R.2). Missing inputs give empty tables.
set -euo pipefail

RUN_DIR="@16_dir1@"

python3 ../integration/integrate.py "$RUN_DIR"
