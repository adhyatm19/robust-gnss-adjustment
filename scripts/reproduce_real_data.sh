#!/usr/bin/env bash
# Reproduce every real-data result in results/rtklib_validation and results/real_network.
#
#   scripts/reproduce_real_data.sh [PYTHON]
#
# Requires git, gcc, make, network access to igs.bkg.bund.de (anonymous HTTPS,
# ~126 MB) and the package installed with the 'rinex' extra:
#   pip install -e ".[rinex,dev]"
# Steps are idempotent: existing downloads and the RTKLIB build are reused.
set -euo pipefail
PY="${1:-python}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export PYTHONDONTWRITEBYTECODE=1

[ -x tools/rtklib/bin/rnx2rtkp-hp ] || scripts/build_rtklib.sh
"$PY" scripts/validate_rtklib_patch.py                                   # Stage 1 checks (GSI pair)
"$PY" scripts/fetch_real_data.py configs/real_data/madrid_2024.yaml      # data + checksums
"$PY" scripts/process_real_network.py configs/real_data/madrid_2024_processing.yaml
"$PY" scripts/build_real_network.py                                      # calibration + networks
for net in main_formal main_calibrated quiet_formal; do                  # Experiment A
    "$PY" -m gnss_adjust adjust --network "results/real_network/networks/$net" \
        --reference "results/real_network/networks/$net/reference_coordinates.csv" \
        --title "Experiment A: Madrid real network ($net)" \
        --output "results/real_network/experiment_A/$net"
done
"$PY" scripts/real_injection_experiments.py --trials 40                  # Experiments B, C
"$PY" scripts/real_antenna_experiment.py                                 # Experiment D
"$PY" scripts/plot_real_results.py                                       # figures
