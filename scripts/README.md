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

⚠️ **Nota de mantenimiento**: Los dos generadores no son interoperables (naming/regex diverge). Ver `docs/` para detalles.

## plot/ — Visualización

```
plot/
├── plot_entropy_paper.jl       # Entropía vs t (FBC/PBC, Plots.jl)
├── plot_entropy_pbc_complete.jl # Merge PBC files, figura Δκ=0.05–0.9 (Plots.jl)
├── plot_entropy_param_sweep.jl  # Entropía genérica multi-param (Plots.jl)
├── plot_ftmle.jl               # λ(t) log-log con power-law (Plots.jl)
├── plot_ensemble_results.jl    # Heatmaps energía modal + termalización (CairoMakie)
└── plot_heatmap_grid.jl        # Grid de heatmaps param×delta (CairoMakie)
```

⚠️ **Nota de backend**: Mitad Plots.jl, mitad CairoMakie. Ver `docs/plotting_backend_diagnostics.md`.

## Notas de desarrollo

- Muchos scripts heredan include-paths de su ubicación anterior en `scripts/`. Si después de reorganizar algún script no corre, revisar líneas con `include("...")` — probablemente necesite un `cd` al directorio anterior o rutas absolutas `joinpath(@__DIR__, ...)`.
- Los scripts de `plot/` que usan `results/figures/...` asumen ejecución desde la raíz del proyecto (`julia --project=. scripts/plot/plot_*.jl`).
- Ver `docs/` para diagnósticos conocidos (duplicación código, HPC workflow breaks, física).
