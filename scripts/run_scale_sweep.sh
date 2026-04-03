#!/bin/bash
# Helper script to run the scale sweep with optional arguments
# Usage:
#   ./scripts/run_scale_sweep.sh              # Default: 10^6, N=32,64,128,256
#   ./scripts/run_scale_sweep.sh 1e7          # TMAX=10^7, default N
#   ./scripts/run_scale_sweep.sh 1e7 "32,64,128,256,512"  # Custom TMAX and N

set -euo pipefail

TMAX="${1:-1e6}"
NS="${2:-32,64,128,256}"

cd "$(dirname "$0")/.."
julia --project=. examples/run_scale_sweep.jl --tmax="$TMAX" --ns="$NS"
