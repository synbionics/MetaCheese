#!/bin/bash
# Step R.2 (final step): interactive HTML report.
# Reads the outputs of all previous steps of both branches and the joint analyses of step R.1
# (missing ones are shown as "not run") and writes <run>/report/report.html, plus
# <run>/report/<run>_report.zip to share it (the report and the pages it links to).
set -euo pipefail

RUN_DIR="@17_dir1@"

python3 ../report/make_report.py "$RUN_DIR"
