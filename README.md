# Diatomic FPUT chain: resonances, thermalization and chaos

Julia code for the numerical study of the **diatomic Fermi–Pasta–Ulam–Tsingou chain**
(alternating disorder in the spring constants, `κ = 1 ± Δκ`, and optionally in the masses,
`m = 1 ± Δm`). The project covers three fronts:

1. **Weakly nonlinear wave theory**: two-branch dispersion relation (acoustic/optical),
   three-wave resonance manifolds ω₋(k₁)+ω₋(k₂)=ω₊(k₃), coupling coefficients Γ, and the
   effective scattering rate R(Δκ).
2. **Dynamics and thermalization**: symplectic integration of trajectories, modal energies,
   spectral entropy S̄(t), thermalization times and interband transfer. Includes
   **random-phase ensembles over a selected band** (see
   [plan_ensamble_fases.md](plan_ensamble_fases.md)) and sweeps over the system size N.
3. **Chaos**: finite-time maximal Lyapunov exponent (ftMLE) via the two-trajectory method
   of Benettin.

The production workflows are meant to run on **HPC with PBS/Torque**; see
[GUIA_HPC_MULTI_N.md](GUIA_HPC_MULTI_N.md).

---

## Repository layout

```
├── src/                    # Core modules
├── scripts/                # Compute and figure pipelines (production)
├── examples/               # Theoretical calculations and exploratory figures
├── configs/
│   ├── production/         # Concrete run configurations (cases)
│   ├── hpc_templates/      # Parameterizable templates for job generation
│   └── tests/              # Smoke tests and quick test configs
├── docs/                   # Supplementary documentation
├── jobs/                   # Generated HPC jobs (.toml versioned, .pbs ignored)
├── results/                # Outputs (mostly git-ignored)
│   ├── data/               # Simulation .jld2 files
│   ├── figures/            # Figures
│   ├── logs/               # Logs and .SUCCESS/.FAILED markers
│   └── scattering_rate/    # R_vs_delta.csv (the only versioned output)
├── supplementary_material/ # Paper PDFs (git-ignored)
├── GUIA_HPC_MULTI_N.md     # Step-by-step guide to the multi-N HPC workflow
├── plan_ensamble_fases.md  # Design of the phase ensemble
├── Project.toml / Manifest.toml
└── README.md
```

---

## Core modules (`src/`)

| File | Module | Contents |
|---|---|---|
| [fput_core.jl](src/fput_core.jl) | `FPUTCore` | `SystemParams`, `make_system` (alternating κ and m), `fput_forces!` (α/β forces with κ scaling, `:fixed` and `:periodic` boundaries), `find_normal_modes` (diagonalization of the dynamical matrix) |
| [fput_fast_runner.jl](src/fput_fast_runner.jl) | `FPUTFastRunner` | `solve_fput`: `SecondOrderODEProblem` integrated with **KahanLi8** (8th-order symplectic) |
| [fput_analysis.jl](src/fput_analysis.jl) | `FPUTAnalysis` | `compute_modal_energies` (mass-weighted projection), `sliding_window_avg` (growing window Δ·t), `spectral_entropy` |
| [config.jl](src/config.jl) | `Config` | `ExperimentConfig` and `PlotConfig` loaded from TOML, plus `default_plot_config` |

The modules are loaded with `include(...)`, not as a registered package.

---

## Compute scripts (`scripts/`)

| Script | What it does |
|---|---|
| [compute_trajectories.jl](scripts/compute_trajectories.jl) | `param_values × delta_values` sweep in time blocks; saves one `.jld2` with modal energies |
| [compute_trajectories_Nsweep.jl](scripts/compute_trajectories_Nsweep.jl) | Same, but sweeping over `N_values` |
| [compute_ensemble.jl](scripts/compute_ensemble.jl) | **Phase ensemble** over a selected band (acoustic/optical): mean and spread of S̄(t), per-branch energies, T_therm. Parallelized over (param, delta) |
| [compute_ensemble_Nsweep.jl](scripts/compute_ensemble_Nsweep.jl) | Ensemble with fixed α and Δκ, sweeping over N (one thread per N) |
| [compute_ftmle.jl](scripts/compute_ftmle.jl) | Benettin ftMLE (two trajectories + renormalization), multiple Δκ in parallel |
| [generate_hpc_jobs.jl](scripts/generate_hpc_jobs.jl) | Generates one `.toml` + `.pbs` per N, with `ppn`/walltime scaled accordingly |
| [check_hpc_status.jl](scripts/check_hpc_status.jl) | Status report for the multi-N jobs (logs, `.SUCCESS`/`.FAILED`, outputs) |

## Figure scripts (`scripts/`)

| Script | Figure |
|---|---|
| [plot_entropy_paper.jl](scripts/plot_entropy_paper.jl) | Official S̄(t) figures for FBC and PBC, with color/style encoding per Δκ |
| [plot_entropy_pbc_complete.jl](scripts/plot_entropy_pbc_complete.jl) | PBC S̄(t) combining the production sweep (Δκ 0.05–0.7) with the Δκ=0.9 run |
| [plot_entropy_param_sweep.jl](scripts/plot_entropy_param_sweep.jl) | Entropy and localization ξ curves for parameter sweeps |
| [plot_ensemble_results.jl](scripts/plot_ensemble_results.jl) | Ensemble diagnostics: modal heatmaps, S̄±σ, global summary, T_th vs Δκ |
| [plot_heatmap_grid.jl](scripts/plot_heatmap_grid.jl) | Grid of modal-energy heatmaps (CairoMakie), configurable via TOML |
| [plot_ftmle.jl](scripts/plot_ftmle.jl) | λ(t) log-log with a fitted power law t^(−δ); supports several overlaid JLD2 files |

## Theory and exploratory scripts (`examples/`)

| Script | Contents |
|---|---|
| [plot_dispersion_relation.jl](examples/plot_dispersion_relation.jl) | ω±(k) for several Δκ |
| [plot_resonance_level_curves.jl](examples/plot_resonance_level_curves.jl) | Level curves of the three-wave resonance residual |
| [coupling_coefficients.jl](examples/coupling_coefficients.jl) | Eigenvectors of the two branches and coupling coefficients |
| [plot_gamma_with_resonance.jl](examples/plot_gamma_with_resonance.jl) | \|Γ\| overlaid on the resonance manifold (provides `compute_gamma`, `omega_branch`) |
| [plot_aao_delta_sweep.jl](examples/plot_aao_delta_sweep.jl) | Δκ sweep of the acoustic+acoustic→optical channel |
| [compute_scattering_rate.jl](examples/compute_scattering_rate.jl) | R(η) as a line integral over the resonance manifold (co-area formula) |
| [plot_entropy_Nsweep.jl](examples/plot_entropy_Nsweep.jl) | Same, with publication styling |
| [plot_thermalization_time.jl](examples/plot_thermalization_time.jl) | T_th vs Δκ with a halo showing the spread across realizations |

---

## Installation

```bash
git clone <repo-url>
cd Codigo
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

`Project.toml` declares: `Colors`, `DifferentialEquations`, `JLD2`, `LaTeXStrings`,
`LinearAlgebra`, `Plots`, `Statistics`.

> **Note**: some scripts use packages that are **not** declared in `Project.toml`.
> Add them to the environment before running those scripts:
>
> ```bash
> julia --project=. -e 'import Pkg; Pkg.add(["CairoMakie","Contour","Interpolations","TOML","Random","Printf","Dates"])'
> ```
>
> - `CairoMakie` — `plot_heatmap_grid.jl`, `plot_ensemble_results.jl` and most of `examples/`
> - `Contour`, `Interpolations` — `examples/compute_scattering_rate.jl`
> - `TOML`, `Random`, `Printf`, `Dates` — stdlibs used by the compute scripts

---

## Usage

### Trajectory simulation

```bash
julia --project=. scripts/compute/compute_trajectories.jl configs/production/periodic_N64.toml
```

### Phase ensemble (production)

```bash
julia --project=. -t 8 scripts/compute/compute_ensemble.jl configs/production/ensemble.toml
julia --project=. scripts/plot/plot_ensemble_results.jl results/data/ensemble_production_100real/ensemble_results_YYYY-MM-DD.jld2
```

### Lyapunov exponent (ftMLE)

```bash
# config.toml, T_max, T_renorm, comma-separated Δκ list
julia --project=. -t 4 scripts/compute/compute_ftmle.jl configs/production/periodic_N64.toml 1e7 200.0 0.1,0.5
julia --project=. scripts/plot/plot_ftmle.jl results/data/ftmle/ftmle_periodic_delta0p1_YYYY-MM-DD.jld2
```

### Paper figures

```bash
julia --project=. scripts/plot/plot_entropy_paper.jl <fbc.jld2> <pbc.jld2>
julia --project=. scripts/plot/plot_entropy_pbc_complete.jl <pbc_0p05-0p7.jld2> <pbc_delta09.jld2>
```

### Multi-N HPC workflow (PBS/TORQUE)

```bash
julia --project=. scripts/hpc/generate_hpc_jobs.jl configs/hpc_templates/ensemble_N_sweep.toml
for job in jobs/ensemble_N*.pbs; do qsub $job; done
julia --project=. scripts/hpc/check_hpc_status.jl
```


---

## Configuration

The simulation TOMLs (`configs/production/`, `configs/hpc_templates/`, `configs/tests/`) have three sections:

```toml
[experiment]
name        = "ensemble_production"
description = "Phase ensemble, Δκ sweep, n_real=10"

[physics]
N              = 64
boundary       = "periodic"        # "periodic" | "fixed"
system_type    = "springs"
nonlinear      = "alpha"           # "alpha" | "beta"
param_values   = [0.1]             # nonlinearity parameter values
delta_values   = [0.1, 0.2, 0.5]   # Δκ sweep
energy_density = 0.445             # E_total = N · energy_density
init_type      = "band_ensemble"   # "mode" | "band_ensemble" (compute_ensemble* only)
branch         = "acoustic"        # branch excited initially
n_real         = 10                # ensemble realizations
seed_base      = 42

[simulation]
TMAX          = 1e6
T_block       = 5e4                # block-wise integration
DT            = 0.05
save_every    = 100
downsample    = 1
entropy_delta = 0.6                # smoothing window for S̄
debug         = true

[output]
base_dir = "results/data/ensemble_production_100real"
```

**Initialization key, per script**: `compute_ensemble.jl` and `compute_ensemble_Nsweep.jl`
read `init_type` (`"mode"` or `"band_ensemble"`), whereas `compute_trajectories.jl` reads
`initial_condition` (default `"low"`) and `compute_trajectories_Nsweep.jl` excites mode
`init_mode`.

**Energy convention**: if the TOML defines `energy_density`, then
`E_total = N · energy_density`; if it only defines `initial_energy`, then
`E_total = initial_energy` (legacy behavior).

Plotting scripts are generally hardcoded; plotting TOMLs are deprecated.

---

## Output format

The compute scripts save a single `.jld2` with two keys: `results` (a vector of
NamedTuples, one per case) and `config` (path of the TOML used). Typical fields of an
ensemble case:

`param`, `Delta`, `scaled_t`, `entropy_mean`, `entropy_std`, `modal_E_mean`,
`E_acoustic_mean`, `E_optical_mean`, `entropy_realizations`, `E_optical_realizations`,
`T_therm_mean`, `T_therm_std`, `T_therm_median`, `T_therm_vec`, `frac_therm`, `n_therm`,
`n_real`, `seed_base`, `branch`, `k_band`, `E_total`.

`compute_ftmle.jl` instead saves flat keys: `t_physical`, `t_cycles`, `lambda`,
`omega_ref`, `delta_k`, `boundary`, `init_mode`, `T_max`, `config_path`.

---

## What is kept out of git

[.gitignore](.gitignore) deliberately excludes the heavy data and figures:

- `results/raw/`, `results/data/`, `results/logs/`, `results/figures/`
- Any `*.jld2`, `*.h5`, `*.jld`
- Any `*.png`, `*.pdf`, `*.svg`, `*.eps` — this includes the PDFs in
  `supplementary_material/`, which exist on disk but are **not** versioned
- `jobs/*.pbs` and `jobs/*.log` (the `.pbs` files are generated; the job `.toml` files are versioned)
- Scheduler output (`*.e*`, `*.o*`), editor/IDE artifacts and `.claude/`

The only versioned output is [results/scattering_rate/R_vs_delta.csv](results/scattering_rate/R_vs_delta.csv).
Reproducing everything else requires re-running the pipelines.

---

## License

MIT
