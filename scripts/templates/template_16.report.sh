#!/bin/bash
# Step R.1 (joint step of both branches): final interactive HTML report.
# Reads the outputs of all previous steps of both branches (missing ones are shown as "not run")
# and writes <run>/report/report.html, a single self-contained file.
set -euo pipefail

RUN_DIR="@16_dir1@"

python3 ../report/make_report.py "$RUN_DIR"
