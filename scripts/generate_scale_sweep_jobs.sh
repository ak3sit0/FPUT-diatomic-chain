#!/bin/bash
set -euo pipefail

# Genera un conjunto de configs y sus jobs PBS para N={32,64,128,256,512}
# 1) Copia el config base
# 2) Ajusta N
# 3) Llama a generate_pbs.jl

BASE_CONFIG="configs/alpha_sweep_N_periodic.toml"
if [[ ! -f "$BASE_CONFIG" ]]; then
  echo "ERROR: no existe $BASE_CONFIG"
  exit 1
fi

# N values para sweep
Ns=(32 64 128 256 512)

for N in "${Ns[@]}"; do
  TARGET_CONFIG="configs/alpha_sweep_N${N}_periodic.toml"
  echo "Generando config $TARGET_CONFIG..."
  cp "$BASE_CONFIG" "$TARGET_CONFIG"
  sed -i "s/^N\s*=.*/N = ${N}/" "$TARGET_CONFIG"
  
  echo "Generando PBS para N=$N..."
  julia --project=. scripts/generate_pbs.jl "$TARGET_CONFIG" \
    --ppn=8 --nodes=1 --mem=64gb --walltime=24:00:00 --queue=workq
  
  echo "-> Generado jobs/$(basename "${TARGET_CONFIG%.toml}").pbs"
  echo
done

echo "Listo. Archivos .toml y .pbs generados para todos los N."
