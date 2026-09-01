# Diatomic FPUT chain: resonances, thermalization and chaos

Julia code for the numerical study of the **diatomic Fermi–Pasta–Ulam–Tsingou chain**
(alternating disorder in the spring constants, `κ = 1 ± Δκ`, and optionally in the masses,
`m = 1 ± Δm`). The project covers three fronts:

1. **Weakly nonlinear wave theory**: two-branch dispersion relation (acoustic/optical),
   three-wave resonance manifolds ω₋(k₁)+ω₋(k₂)=ω₊(k₃), coupling coefficients Γ, and the
   effective scattering rate R(Δκ).
2. **Dynamics and thermalization**: symplectic integration of trajectories, modal energies,
   spectral entropy S̄(t), thermalization times and interband transfer. Includes
   **random-phase ensembles over a selected band** and sweeps over the system size N.
3. **Chaos**: finite-time maximal Lyapunov exponent (ftMLE) via the two-trajectory method
   of Benettin.

The production workflows are meant to run on **HPC with PBS/Torque**.

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
├── test/                   # Unit tests (julia --project=. test/runtests.jl)
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
| [config.jl](src/config.jl) | `Config` | `ExperimentConfig` and `PlotConfig` loaded from TOML, `default_plot_config`, and `resolve_config_path`/`as_experiment_config` for results files that recorded a since-moved config |
| [experiment.jl](src/experiment.jl) | `Experiment` | Reads a sweep TOML into one `SweepSpec`. Stores the *policy* (ε vs E_total, `scaled_t_max` vs `TMAX`, `n_samples` vs `save_every`); `resolve_budget` turns it into concrete numbers per N |
| [case_setup.jl](src/case_setup.jl) | `CaseSetup` | `build_case` → a fully specified initial value problem (system, normal modes, bands, energy budget); `initial_condition` draws the starting state |
| [sweep_driver.jl](src/sweep_driver.jl) | `SweepDriver` | Runs a task list, survives individual failures, saves the surviving results. Parallelism policy is an argument, not a hardcoded choice |
| [plotting_utils.jl](src/plotting_utils.jl) | `PlottingUtils` | Single source for palettes/linestyles (`PALETTE_DELTA`, `BLUES_STOPS`, `SPEC_*`), the `Curve` type that figure scripts build, and shared helpers (`entropy_series`, `prepare_ts`, `logdownsample`) |
| [plot_style.jl](src/plot_style.jl) | `PlotStyle` | The look of every Plots.jl figure: `apply_style!`, `logplot`, `draw!`, `save_fig` |

The modules are loaded with `include(...)`, not as a registered package.

---

## Compute scripts (`scripts/`)

Each one takes a single TOML and sweeps **every combination of `N × param_values × delta_values`**
found in it. A sweep over system size is not a separate script: just list more than one
value in `N_values`.

| Script | What it does |
|---|---|
| [compute_trajectories.jl](scripts/compute/compute_trajectories.jl) | One deterministic trajectory per combination, all the energy in a single mode. Saves modal energies and S̄(t) |
| [compute_ensemble.jl](scripts/compute/compute_ensemble.jl) | **Phase ensemble** over a selected band (acoustic/optical): mean and spread of S̄(t), per-branch energies, thermalization times. `n_real` sets the number of realizations |
| [compute_ftmle.jl](scripts/compute/compute_ftmle.jl) | Benettin ftMLE (two trajectories + renormalization), multiple Δκ in parallel |

## Figure scripts (`scripts/`)

| Script | Figure |
|---|---|
| [plot_entropy_param_sweep.jl](scripts/plot/plot_entropy_param_sweep.jl) | **The base figure**: S̄(t) on log time, one curve per Δκ. Every other time-series figure is a variation of this one |
| [plot_entropy_paper.jl](scripts/plot/plot_entropy_paper.jl) | Same, with hand-picked colors and line widths for the paper. FBC and PBC; extra PBC files merge into a combined Δκ 0.05–0.9 figure |
| [plot_entropy_size_sweep.jl](scripts/plot/plot_entropy_size_sweep.jl) | Same, but one curve per **N** at fixed (param, Δκ), taking one run from each of several files |
| [plot_ftmle.jl](scripts/plot/plot_ftmle.jl) | λ(t) log-log with a fitted power law t^(−δ); supports several overlaid JLD2 files |
| [plot_heatmap_grid.jl](scripts/plot/plot_heatmap_grid.jl) | Grid of modal-energy heatmaps (CairoMakie), configurable via TOML |

## Theory and exploratory scripts (`examples/`)

| Script | Contents |
|---|---|
| [plot_dispersion_relation.jl](examples/dispersion/plot_dispersion_relation.jl) | ω±(k) for several Δκ |
| [plot_resonance_level_curves.jl](examples/resonance/plot_resonance_level_curves.jl) | Level curves of the three-wave resonance residual |
| [coupling_coefficients.jl](examples/coupling/coupling_coefficients.jl) | 3×3 grid of \|Γ_σ1σ2σ3\|; `resonance=true` overlays the resonance manifold |
| [plot_aao_delta_sweep.jl](examples/coupling/plot_aao_delta_sweep.jl) | Δκ sweep of the acoustic+acoustic→optical channel |
| [compute_scattering_rate.jl](examples/coupling/compute_scattering_rate.jl) | R(η) as a line integral over the resonance manifold (co-area formula) |

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
> - `CairoMakie` — `plot_heatmap_grid.jl` and most of `examples/`
> - `Contour`, `Interpolations` — `examples/coupling/compute_scattering_rate.jl`
> - `TOML`, `Random`, `Printf`, `Dates` — stdlibs used by the compute scripts

---

## Usage

Every run is **two commands**: one to compute, one to plot. The compute step reads a TOML
and writes a `.jld2`; the plot step reads that `.jld2` and writes a PDF and a PNG. You never
have to edit Julia code to change a sweep — change the TOML.

### 1. Trajectories → entropy figure

```bash
julia --project=. scripts/compute/compute_trajectories.jl configs/production/periodic_N64.toml
julia --project=. scripts/plot/plot_entropy_param_sweep.jl results/data/sweep_results_YYYY-MM-DD.jld2
```

To sweep system sizes instead, put several values in `N_values` in the TOML — same script,
same command. The figure then labels each curve with its N.

### 2. Phase ensemble

```bash
julia --project=. -t 8 scripts/compute/compute_ensemble.jl configs/production/ensemble.toml
```

`-t 8` gives Julia 8 threads. Where they are used is decided for you: with `n_real > 1` the
realizations run in parallel inside each case; with `n_real = 1` the cases themselves run in
parallel. Either way the result is reproducible bit for bit from `seed_base`.

### 3. Lyapunov exponent (ftMLE)

```bash
# arguments: config.toml, T_max, T_renorm, comma-separated Δκ list
julia --project=. -t 4 scripts/compute/compute_ftmle.jl configs/production/periodic_N64.toml 1e7 200.0 0.1,0.5
julia --project=. scripts/plot/plot_ftmle.jl results/data/ftmle/ftmle_periodic_delta0p1_YYYY-MM-DD.jld2
```

### 4. Comparing system sizes at one Δκ

```bash
julia --project=. scripts/plot/plot_entropy_size_sweep.jl --param 0.1 --delta 0.1 <N64.jld2> <N128.jld2> <N256.jld2>
```

### 5. Paper figures

```bash
# third and later arguments add the merged Δκ 0.05–0.9 figure
julia --project=. scripts/plot/plot_entropy_paper.jl <fbc.jld2> <pbc.jld2> [<pbc_delta09.jld2> ...]
```

### Checking your setup

```bash
julia --project=. test/runtests.jl
```

Should report all tests passing. It checks the physics (energy conservation, FPUT recurrence,
dispersion) and the config → simulation pipeline, and takes well under a minute.

### On HPC (PBS/TORQUE)

Write one `.pbs` per job that calls the same commands as above, and submit them:

```bash
for job in jobs/ensemble_N*.pbs; do qsub $job; done
```

The job generators used for the multi-N campaign now live in `deprecated/hpc/` and are
not versioned.


---

## Configuration

The simulation TOMLs (`configs/production/`, `configs/hpc_templates/`, `configs/tests/`) have three sections:

```toml
[experiment]
name        = "ensemble_production"
description = "Phase ensemble, Δκ sweep, n_real=10"

[physics]
N_values       = [64]              # list several values to sweep system size
boundary       = "periodic"        # "periodic" | "fixed"
system_type    = "springs"         # "springs" (Δκ) | "masses" (Δm)
nonlinear      = "alpha"           # "alpha" | "beta"
param_values   = [0.1]             # nonlinearity parameter values
delta_values   = [0.1, 0.2, 0.5]   # Δκ sweep
energy_density = 0.00695           # E_total = N · energy_density
init_type      = "band_ensemble"   # "mode" | "band_ensemble"
branch         = "acoustic"        # band excited initially (band_ensemble only)
n_real         = 10                # ensemble realizations
seed_base      = 42

[simulation]
scaled_t_max  = 6e4                # observation window in cycles (TMAX derived per N)
n_blocks      = 20                 # block-wise integration
DT            = 0.05
n_samples     = 2000               # stored samples, held constant across N
downsample    = 1
entropy_delta = 0.6                # smoothing window for S̄
debug         = false

[output]
base_dir = "results/data/ensemble_production_100real"
```

All three compute scripts read the **same keys with the same meaning**. Where a key can be
specified two ways, either is valid — pick one:

| Choice | Options | What each means |
|---|---|---|
| Energy | `energy_density` **or** `initial_energy` | ε fixed, so `E_total = ε·N` (recommended for N sweeps) **or** `E_total` fixed, so ε falls as 1/N |
| Time budget | `scaled_t_max` **or** `TMAX` | Observation window in cycles, so `TMAX` grows with N **or** a fixed physical time, so large N see less of the dynamics |
| Sampling | `n_samples` **or** `save_every` | Fixed number of stored samples, keeping memory ∝ N **or** a fixed stride, making memory ∝ N² |
| Blocks | `n_blocks` **or** `T_block` | Split `TMAX` into this many blocks **or** use blocks of this length |

**Which mode gets excited** (`init_type = "mode"`): by default the lowest physical mode —
mode 1 for `fixed`, mode 2 for `periodic`, because under periodic boundaries mode 1 is the
uniform translation and does not oscillate. Set `init_mode` to override; the run reports it
when you do. Under `init_type = "band_ensemble"` a whole band is excited instead, chosen by
`branch` (or by `k_band_start`/`k_band_end` for a manual range).

Plotting scripts are generally hardcoded; plotting TOMLs are deprecated.

---

## Output format

The compute scripts save a single `.jld2` with two keys: `results` (a vector of
NamedTuples, one per case) and `config` (path of the TOML used). Every case carries
`N`, `param`, `Delta`, `scaled_t`, `E_total`, `energy_density`, `omega_ref`, `TMAX` and
`save_every`, so a result is self-describing without going back to the TOML.

A trajectory case adds `modal_E`, `entropy`, `t_abs`. An ensemble case adds
`entropy_mean`, `entropy_std`, `modal_E_mean`, `E_acoustic_mean`, `E_optical_mean`,
`entropy_realizations`, `E_optical_realizations`, `T_therm_mean`, `T_therm_std`,
`T_therm_median`, `T_therm_vec`, `frac_therm`, `n_therm`, `n_real`, `seed_base`,
`branch`, `k_band`.

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

The `jobs/` directory is local-only: `.gitkeep` marks the directory, but generated
`.pbs` and `.toml` files (written by the HPC job generators) are never versioned. They
serve as an audit trail on your machine for which configurations were submitted; to clean
up, run `git clean -fd jobs/` (after backing up any `.pbs` scripts you want to keep).

---

## License

MIT
