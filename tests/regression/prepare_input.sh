#!/usr/bin/env bash
# Input of the regression test: samples 1, 2 and 3 of each group, first 5M read
# pairs of each, written into input/test_regression/.
#
# Usage (repository root, on the host):
#   bash tests/regression/prepare_input.sh                  # stream from ENA
#   bash tests/regression/prepare_input.sh <FULL_DATA_DIR>  # local full dataset
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

export OUT_DIR="$REPO_ROOT/input/test_regression"
export READ_PAIRS=5000000
export SAMPLES="\
D_GP_01 SRR31729974 SAMN45828309 cheese01 PDO
D_GP_02 SRR31729973 SAMN45828310 cheese02 PDO
D_GP_03 SRR31729954 SAMN45828336 cheese03 PDO
ND_GP_1 SRR31729966 SAMN45828325 Ncheese1 non-PDO
ND_GP_2 SRR31729959 SAMN45828331 Ncheese2 non-PDO
ND_GP_3 SRR31729958 SAMN45828332 Ncheese3 non-PDO"

exec bash "$REPO_ROOT/examples/small_PRJNA1198075/prepare_input.sh" "$@"
