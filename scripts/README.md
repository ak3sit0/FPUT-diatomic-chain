# Scripts — Simulación FPUT y análisis

Organización temática: **compute/** (simulaciones), **hpc/** (generación de jobs y monitoreo), **plot/** (visualización).

## compute/ — Trayectorias y análisis dinámico

```
compute/
├── calibrate_cost.jl           # Benchmark coste máquina C (ns/step/sitio)
├── compute_trajectories.jl     # Barrido param×delta, energía modal por trayectoria
├── compute_trajectories_Nsweep.jl # Igual, con barrido adicional en N
├── compute_ensemble.jl         # Ensamble n_real realizaciones + T_therm
├── compute_ensemble_Nsweep.jl  # Ensamble con barrido N (per-delta)
└── compute_ftmle.jl            # Exponente Lyapunov finito (ftMLE, Benettin)
```

## hpc/ — Cluster workflow

```
hpc/
├── generate_pbs_nsweep.jl      # PBS jobs para compute_ensemble_Nsweep.jl
├── generate_hpc_jobs.jl        # PBS jobs para compute_ensemble.jl
└── check_hpc_status.jl         # Monitorea logs y reporta SUCCESS/FAILED
```

Workflow: `generate_*.jl` → `qsub` → `check_hpc_status.jl` → `compute_ensemble*.jl` output.

`check_hpc_status.jl` descubre jobs directamente desde `jobs/*.pbs` (no asume naming ni N_values fijo), así que funciona con ambos generadores indistintamente.

## plot/ — Visualización

```
plot/
├── plot_entropy_paper.jl       # Entropía vs t: FBC, PBC y PBC-completo (Plots.jl)
├── plot_entropy_param_sweep.jl  # Entropía genérica multi-param (Plots.jl)
├── plot_ftmle.jl               # λ(t) log-log con power-law (Plots.jl)
├── plot_ensemble_results.jl    # Heatmaps energía modal + termalización (CairoMakie)
└── plot_heatmap_grid.jl        # Grid de heatmaps param×delta (CairoMakie)
```

⚠️ **Nota de backend**: Mitad Plots.jl, mitad CairoMakie. Ver `docs/plotting_backend_diagnostics.md`.

## Notas de desarrollo

- Los scripts de `plot/` que usan `results/figures/...` asumen ejecución desde la raíz del proyecto (`julia --project=. scripts/plot/plot_*.jl`).
- Paletas, estilos y helpers de entropía (`entropy_series`, `prepare_ts`, `logdownsample`,
  `load_results_by_delta`, `cyc`) viven en `PlottingUtils` (`src/plotting_utils.jl`). Usa
  siempre `cyc(palette, i)` en vez de `palette[i]` para indexar paletas.
- `plot_entropy_pbc_complete.jl` se fusionó en `plot_entropy_paper.jl`: era el mismo script con
  otra spec de Δκ y un loader multi-archivo. Pasar un tercer argumento genera la figura extendida.
- `derive_band_indices`/`band_phase_ic` viven en `FPUTCore` (`src/fput_core.jl`); el loop de integración por bloques (usado en los 4 scripts `compute_*.jl`) vive en `BlockIntegration.integrate_in_blocks` (`src/block_integration.jl`).
- Ver `docs/` para diagnósticos conocidos (duplicación código restante, backend de plotting, física).
