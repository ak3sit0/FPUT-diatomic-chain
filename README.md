# FPUT Energy Localization & Entropy Analysis (Literate & Functional)

Simulador de cadenas de Fermi-Pasta-Ulam-Tsingou (FPUT) rediseñado bajo principios de **SICP**, **Composing Programs** y **Literate Programming**.

## Arquitectura (Layers)
El repositorio se organiza en capas lógicas para minimizar el estado mutable y maximizar la legibilidad:

1.  **FPUTCore (`src/fput_core.jl`)**: Definiciones físicas puras (masas, resortes, fuerzas Hamiltonianas).
2.  **FPUTFastRunner (`src/fput_fast_runner.jl`)**: Integración numérica funcional usando solvers de alta precisión (DifferentialEquations.jl).
3.  **FPUTAnalysis (`src/fput_analysis.jl`)**: Algoritmos de análisis (Transformada modal, Entropía Espectral).
4.  **Orquestadores (`examples/`)**: Scripts literarios que cuentan la historia de la simulación y visualización.

## Quick Start

### Installation

```bash
git clone <repo-url>
cd Codigo
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

### Running a Simulation

```bash
julia --project=. examples/fput_literate_sim.jl configs/templates/test_quick.toml
```

### Generating Plots

```bash
julia --project=. examples/fput_literate_viz.jl results/raw/refactored_results_YYYY-MM-DD.jld2
```

### HPC Workflow (offline job generation)

Generate `.pbs` job scripts without submitting them. This writes paired `jobs/*.pbs` and `jobs/*.toml` files so you can inspect and submit manually:

```bash
julia --project=. scripts/generate_hpc_jobs.jl --template configs/templates/test_quick.toml --sweep N
```

Then submit manually (example):

```bash
for job in jobs/*.pbs; do qsub "$job"; done
```

Notes:
- The generator mutates a copy of the template TOML and writes a PBS wrapper that calls `examples/fput_literate_sim.jl` with the job-specific TOML.
- Simulation outputs are saved under the `base_dir` defined in each TOML. The saved `.jld2` file contains a primitive `config` Dict (keys: `N`, `DeltaK`, `DeltaM`, `alpha`, `beta`, `boundary`) to avoid type-reconstruction warnings when loading on different machines.
- Use `examples/fput_literate_viz.jl` to produce separate PDF outputs (`_heatmap.pdf`, `_entropy.pdf`, `_resonance.pdf`).

## Estructura del Proyecto
- `src/`: Motores de cálculo (FPUTCore, FastRunner, Analysis).
- `configs/`: Plantillas TOML descriptivas.
- `examples/`: Scripts literarios principales (Simulación y Visualización).
- `results/`: Salidas `.jld2` y reportes visuales en PDF.
- `old_repo_archive.tar.gz`: Resguardo de la versión anterior (tesis original).

## Filosofía
Este código busca la **elegancia y simplicidad**. Se han eliminado docenas de archivos redundantes y dependencias pesadas en favor de funciones puras y despachos por tipos claros en Julia.


[output]
entropy_dir  = "results/figures/entropy"
xi_dir       = "results/figures/xi"
heatmaps_dir = "results/figures/heatmaps"
EOF

# Generate entropy & localization plots
julia --project=. examples/plot_entropy_data.jl configs/plot_my_experiment.toml


Unified plotting
----------------
This repository now provides a unified plotting workflow driven by a single
TOML template: `configs/templates/master_plot.toml`. Use `examples/plot_dynamics.jl`
for entropy/ξ/time-series plots, `examples/plot_heatmap_grid.jl` for heatmaps and
`examples/plot_resonance_level_curves.jl` for resonance visualizations. A small
test runner is available at `examples/run_plot_compat_tests.jl` to smoke-test
these scripts against the existing TOML configs. Example:

```bash
julia --project=. examples/run_plot_compat_tests.jl
```
# Generate energy heatmap grid
julia --project=. examples/plot_heatmap_grid.jl configs/plot_my_experiment.toml
```

## Directory Structure

```
├── configs/                   # TOML configuration files (experiments & plots)
├── examples/                  # End-to-end Julia scripts
│   ├── compute_modal_energies.jl    (main simulation runner)
│   ├── plot_entropy_data.jl         (entropy & localization analysis)
│   ├── plot_heatmap_grid.jl         (energy heatmaps)
│   └── make_animation.jl            (FPUT chain visualization)
├── src/                       # Library modules
│   ├── config.jl              (ExperimentConfig & PlotConfig structs)
│   ├── parameters.jl          (physical constants, make_k_m)
│   ├── simulation_runner.jl   (physics integrator wrapper)
│   ├── dynamical_matrix.jl    (normal mode decomposition)
│   ├── fput_equations.jl      (Hamiltonian formulation)
│   ├── integrator.jl          (Störmer–Verlet scheme)
│   ├── energy_analysis.jl     (entropy & mode energy calculations)
├── scripts/                   # Utilities
│   └── generate_pbs.jl        (auto-generate HPC job files)
├── results/                   # Output directory (not in git)
│   ├── raw/                   (simulation data, .jld2 files)
│   ├── figures/               (plots)
│   └── logs/                  (job logs)
├── tests/                     # Unit tests
├── Project.toml               # Julia package manifest
└── Manifest.toml
```

## Configuration Format

TOML files contain three sections:


nonlinear         = "alpha"               # or "beta"
param_values      = [0.1, 0.2]            # α or β values
delta_values      = [0.1, 0.3, 0.5]       # ΔK, ΔM values
initial_condition = "low"                 # or "high"
initial_energy    = 0.45

[simulation]
TMAX       = 1e7
T_block    = 1e5       # Checkpoint interval
DT         = 0.05      # Time step
save_every = 1000
downsample = 2
debug      = false

[output]
base_dir = "results/raw"
```

## Key Modules

| Module | Purpose |
|--------|---------|
| `Config` | Load/validate TOML files; save results with reproducibility metadata |
| `Parameters` | Physical constants (DT) and spring-mass builder |
| `SimulationRunner` | FPUT ODE integrator (Störmer–Verlet) |
| `DynamicalMatrix` | Normal modes & frequencies |
| `EnergyAnalysis` | Modal energy, spectral entropy, localization (ξ) |

## Results & Metadata

Every `.jld2` output file contains:
- Simulation data (energies, time series)
- **Metadata** (TOML config, git commit hash, timestamp, hostname, Julia version)

This enables **full reproducibility**: rerun with the same config → identical results.

## Testing

```bash
julia --project=. tests/test_checkpoint.jl
```

## License

MIT

## References

- Localization parameter (ξ): *Flach & Gorning (2004)*
- FPUT dynamics: *Fermi, Pasta, Ulam, Tsingou*
