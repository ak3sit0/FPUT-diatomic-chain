"""
    compute_trajectories.jl

Deterministic trajectory sweep over `N × param × Δ`: one run per combination,
all the energy in a single mode.

Replaces `compute_trajectories.jl` + `compute_trajectories_Nsweep.jl`. The
N sweep was never a different computation — it is this script with more than one
entry in `N_values`.

Usage:
  julia --project=. scripts/compute/compute_trajectories.jl <config.toml>
"""

using JLD2, Statistics, LinearAlgebra
include("../../src/fput_core.jl");        using .FPUTCore
include("../../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../../src/fput_analysis.jl");    using .FPUTAnalysis
include("../../src/block_integration.jl");using .BlockIntegration
include("../../src/experiment.jl");       using .Experiment
include("../../src/case_setup.jl");       using .CaseSetup
include("../../src/sweep_driver.jl");     using .SweepDriver

const USAGE = "Usage: julia --project=. scripts/compute/compute_trajectories.jl <config.toml>"

"""Integrate one case and return its result entry."""
function run_case(spec, task)
    label = "N=$(task.N) p=$(task.param) Δ=$(task.delta)"
    println("[$label] ε=$(round(spec_epsilon(spec, task.N); sigdigits=4))")

    case   = build_case(spec, task.N, task.param, task.delta)
    q0, v0 = initial_condition(case, spec, spec.seed_base)
    b      = case.budget

    result = integrate_in_blocks(case.sp, q0, v0, case.N, case.freq, case.V, case.m,
                                 case.ref_mode;
                                 TMAX = b.TMAX, T_block = b.T_block, DT = spec.DT,
                                 save_every = b.save_every, downsample = spec.downsample,
                                 track_abs_time = true, debug = spec.debug, label = label)

    entropy = FPUTAnalysis.spectral_entropy(result.modal_E, spec.entropy_delta)
    println("[$label] done; nt=$(length(result.scaled_t))  S_final=$(round(entropy[end]; digits=4))")

    (; N = case.N, param = case.param, Delta = case.delta,
       scaled_t = result.scaled_t, t_abs = result.t_abs,
       modal_E = result.modal_E, entropy = entropy,
       E_total = b.E_total, energy_density = b.E_total / case.N,
       omega_ref = case.omega_ref, TMAX = b.TMAX, save_every = b.save_every,
       init_type = spec.init_type, excited = case.excited)
end

spec_epsilon(spec, N) = isnothing(spec.energy_density) ? spec.E_total / N : spec.energy_density

function main()
    isempty(ARGS) && error(USAGE)
    spec  = parse_spec(ARGS[1])
    tasks = sweep_tasks(spec)

    println("Sweep: $(length(tasks)) cases = $(length(spec.N_values)) N × " *
            "$(length(spec.param_values)) param × $(length(spec.delta_values)) Δ\n")

    run_sweep(t -> run_case(spec, t), tasks;
              outfile     = output_path(spec, "sweep_results"),
              config_path = ARGS[1],
              parallel    = :threads)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
