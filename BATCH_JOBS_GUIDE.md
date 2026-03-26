# 🚀 Batch Job Submission Guide

Comprehensive FPUT parameter sweeps organized in 3 batch PBS jobs.

## Overview

**7 simulation configurations** across spring systems with different:
- Boundary conditions (fixed / periodic)
- Initial modes (acoustic / optical)  
- Nonlinearity type (α / β)

**Organized into 3 PBS jobs** to minimize cluster queue load:

| Job | Configs | Expected Duration |
|-----|---------|------------------|
| `fixed_bc` | 4 (acoustic/optical × α/β) | ~200-300 hours (16 PPNs) |
| `periodic_acoustic` | 2 (α/β) | ~100-150 hours (12 PPNs) |
| `periodic_optical_beta` | 1 (β only, α done separately) | ~50-100 hours (8 PPNs) |

---

## Quick Launch

```bash
# Submit all 3 jobs to cluster
qsub jobs/fixed_bc.pbs
qsub jobs/periodic_acoustic.pbs
qsub jobs/periodic_optical_beta.pbs

# Monitor
qstat | grep fixed_bc
qstat | grep periodic

# Check logs (ongoing)
tail -f results/logs/*.log
```

---

## Configurations Details

### Job 1: Fixed Boundary Conditions

**Runs:**
1. `fixed_acoustic_alpha.toml` → mode 1, α=0.1
2. `fixed_acoustic_beta.toml` → mode 1, β=0.1
3. `fixed_optical_alpha.toml` → mode 33, α=0.1
4. `fixed_optical_beta.toml` → mode 33, β=0.1

**Output:**
```
results/raw/
  fixed_acoustic_alpha/
  fixed_acoustic_beta/
  fixed_optical_alpha/
  fixed_optical_beta/
```

### Job 2: Periodic BC, Acoustic Modes

**Runs:**
1. `periodic_acoustic_alpha.toml` → mode 2, α=0.1
2. `periodic_acoustic_beta.toml` → mode 2, β=0.1

**Output:**
```
results/raw/
  periodic_acoustic_alpha/
  periodic_acoustic_beta/
```

### Job 3: Periodic BC, Optical Mode (β only)

**Runs:**
1. `periodic_optical_beta.toml` → mode 33, β=0.1

**Output:**
```
results/raw/
  periodic_optical_beta/
```

**Note:** `periodic_optical_alpha` (α+optical) done separately in `alpha_optical_N64_periodic` (already submitted).

---

## Common Parameters (All Configs)

```
N = 64
ΔK sweep = [0.0, 0.01, 0.1, 0.12, 0.15, 0.17, 0.2, 0.3, 0.6, 0.8, 0.9]
E₀ = 0.45
TMAX = 10⁸
T_block = 10⁶
DT = 0.05
System type = springs
```

---

## Regenerating PBS Files

If you modify a group config, regenerate its PBS:

```bash
# Modify and save jobs/group_fixed_bc.toml, then:
julia scripts/batch_jobs.jl jobs/group_fixed_bc.toml --ppn=16 --walltime=72:00:00

# Or with custom parameters:
julia scripts/batch_jobs.jl jobs/group_fixed_bc.toml \
    --ppn=20 \
    --mem=256gb \
    --walltime=120:00:00 \
    --queue=gpu
```

---

## Output Organization Flow

After all jobs complete:

```bash
# 1. Verify all results saved
ls results/raw/*/results_springs_*.jld2

# 2. Create plot configurations
# Example: plot only specific ΔK values
cat > configs/plot_acoustic_comparison.toml << 'EOF'
[data]
input_file = "results/raw/fixed_acoustic_alpha/results_*.jld2"

[filter]
delta_values = [0.0, 0.3, 0.6, 0.9]  # subset of ΔK

[output]
entropy_dir = "results/figures/entropy/acoustic_fixed"
EOF

# 3. Run plotting
julia --project=. examples/plot_entropy_data.jl configs/plot_acoustic_comparison.toml
```

---

## Troubleshooting

**Q: One config fails partway through?**
- PBS logs are in `results/logs/`  
- Checkpoints are in `results/raw/*/checkpoints/`
- Re-submit just that job config if needed

**Q: Need to run with different parameters?**
- Edit the TOML (`configs/fixed_acoustic_alpha.toml` etc.)
- Regenerate PBS: `julia scripts/batch_jobs.jl jobs/group_NAME.toml`
- Resubmit: `qsub jobs/NAME.pbs`

**Q: How long will this take?**
- Depends on cluster load + node performance
- Expect 10-14 days total for all 3 jobs (sequential)
- Can submit all 3 simultaneously to run in parallel (~5-7 days)

---

## Advanced: Custom Batch Job

Create a new grouping:

```toml
# jobs/group_custom.toml
[job]
name = "my_custom_sweep"
description = "My specific experiment"

[[simulations]]
config = "configs/fixed_acoustic_alpha.toml"

[[simulations]]
config = "configs/periodic_optical_beta.toml"

[[simulations]]
config = "configs/my_other_config.toml"
```

Then:
```bash
julia scripts/batch_jobs.jl jobs/group_custom.toml --ppn=16
qsub jobs/my_custom_sweep.pbs
```
