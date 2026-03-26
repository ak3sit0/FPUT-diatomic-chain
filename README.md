# FPUT Energy Localization & Entropy Analysis

A Julia toolkit for simulating Fermi-Pasta-Ulam-Tsingou (FPUT) chains and analyzing modal energy dynamics through spectral entropy and localization parameters.

## Features

- **Configurable FPUT simulations** via TOML files (N, boundary conditions, nonlinearity type, energy)
- **Automatic PBS job generation** for HPC cluster submission
- **Post-processing pipeline** with entropy & localization analysis
- **Reproducible workflows** with embedded metadata (git hash, timestamp, config file)
- **Publication-ready visualizations** (heatmaps, entropy plots, animations)

## Quick Start

### Installation

```bash
git clone <repo-url>
cd Codigo
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

### Running a Simulation

1. **Create a configuration file** (copy & modify an existing one):
   ```bash
   cp configs/alpha_sweep_periodic.toml configs/my_experiment.toml
   # Edit TOML: change N, param_values, delta_values, TMAX, etc.
   ```

2. **Run locally**:
   ```bash
   julia --project=. --threads=8 examples/compute_modal_energies.jl configs/my_experiment.toml
   # Results save to: results/raw/my_experiment/
   ```Caso 1: Excito acústico → La energía se queda atrapada (Sticky States / Bloqueo por Gap).

3. **Or submit to cluster**:
   ```bash
   julia scripts/generate_pbs.jl configs/my_experiment.toml --ppn=16 --walltime=48:00:00
   qsub jobs/my_experiment.pbs
   ```

### Generating Plots

```bash
# Create a plot config
cat > configs/plot_my_experiment.toml << 'EOF'
[data]
input_file = "results/raw/my_experiment/results_springs_*.jld2"

[output]
entropy_dir  = "results/figures/entropy"
xi_dir       = "results/figures/xi"
heatmaps_dir = "results/figures/heatmaps"
EOF

# Generate entropy & localization plots
julia --project=. examples/plot_entropy_data.jl configs/plot_my_experiment.toml

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

```toml
[experiment]
name        = "my_sim"
description = "Optional description"

[physics]
N                 = 32                    # System size
boundary          = "periodic"            # or "fixed"
system_type       = "springs"
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
