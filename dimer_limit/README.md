# Δκ → 1: the dimer limit

Isolated workspace for the diatomic FPUT chain (`κ = 1 ± Δκ`) as Δκ approaches, and reaches, 1.
At Δκ = 1 the soft springs vanish (κ₂ = 0) and the chain is N/2 **independent dimers**.

It reuses `../src` and `../scripts` unchanged; only what is specific to this limit lives here.

```
dimer_limit/
├── spectrum.jl      preflight: bands, clock, time budget, distance to the cubic turnover
├── exact_dimers.jl  exact Δκ = 1 reference (decoupled dimers), checked against the full chain
├── configs/         two sweeps for the pipeline
└── data/            outputs (created on demand, git-ignored)
```

## What changes as Δκ → 1

Verified symbolically (SymPy) and against `find_normal_modes`; `spectrum.jl` reproduces the table.

| Quantity | Closed form | Consequence |
|---|---|---|
| Acoustic band edge | ω_ac,max = √(2(1−Δκ)) | band collapses to a point |
| Optical band edge | ω_op,min = √(2(1+Δκ)); ω_op(0) = 2 | band stays at ≈ 2 |
| Acoustic clock (PBC) | N·ω₂ → 2π√(1−Δκ²) | ω₂ is 22× smaller at Δκ=0.999 than at 0.1, so a fixed `scaled_t_max` costs 22× more |
| Δκ = 1 | flat bands at ω = 0 and 2 | 32 degenerate zero modes; no clock, modal basis arbitrary |

## Usage (from the repository root)

```bash
julia --project=. dimer_limit/spectrum.jl                # preflight table, ~seconds
julia --project=. dimer_limit/exact_dimers.jl            # Δκ = 1 control, prints PASS/FAIL

julia --project=. -t 4 scripts/compute/compute_trajectories.jl dimer_limit/configs/periodic_N64_prod_energy.toml
julia --project=. -t 4 scripts/compute/compute_trajectories.jl dimer_limit/configs/periodic_N64_dk0.999.toml

julia --project=. scripts/plot/plot_entropy_param_sweep.jl dimer_limit/data/prod_energy/sweep_results_<date>.jld2
```

`spectrum.jl` takes optional arguments: `N alpha E_total scaled_t_max target_strain`.
Figures land in `results/figures/entropy/`, not here (the plot scripts choose their own output).

## Decisions worth knowing

- **Δκ = 1 is not a pipeline case.** `compute_trajectories.jl` refuses it on purpose (`ref_frequency`,
  ω ≈ 5e-8): with no acoustic band there is no `scaled_t`. `exact_dimers.jl` handles it in dimer
  coordinates and confirms that every dimer energy is conserved, so T_th = ∞ there. That is the endpoint
  the Δκ < 1 runs approach.
- **The energy has to follow Δκ.** Mode 2 moves the dimers rigidly, so it stores its energy in the soft
  bonds and the initial strain grows as 1/√(1−Δκ) at fixed E. At the production energy (E_total = 0.445),
  the largest α·|d| at t = 0 is 0.07 (Δκ=0.9), 0.23 (0.99) and 0.74 (0.999); the cubic potential turns
  over at α·d = 1. Hence two configs: the production energy for Δκ ∈ {0.9, 0.99}, and E_total = 0.044 for
  0.999 so that it sits at the same distance from the turnover as 0.99. The last column of `spectrum.jl`
  gives the E that fixes any target strain if you prefer a fixed-strain protocol over fixed ε.
- **Single-trajectory path, PBC only.** Mode 2 at the production energy (`periodic_N64.toml`), not the
  random-phase ensemble, which uses a different energy density. Under fixed BC the Δκ = 1 point differs
  (31 zero modes plus two edge oscillators at ω = √2), so it is not covered.
- **Separate `base_dir` per config**, because output files are named by date only and a second run the
  same day into the same directory would overwrite the first.

## Caveat, not analyzed

In a throwaway test, `exact_dimers.jl` with random dimer velocities (σ = 0.1) and Δκ = 0.999 instead of 1
diverged (centre-of-mass displacement 3.8·10³ against a ballistic bound of about 2·10²). Consistent with
the soft bonds stretching until the cubic potential turns over, but this was not checked. Worth keeping in
mind before exciting anything but a low-velocity-dispersion state near the limit.
