"""
scripts/compute_trajectories_Nsweep.jl

Barrido sobre N_values definido en un TOML. Guarda un JLD2 con entradas
por cada N: (N, scaled_t, modal_E, entropy).

Usage:
  julia --project=. scripts/compute_trajectories_Nsweep.jl configs/tests/quick_smoke_test.toml
"""

using TOML, JLD2, Dates, Statistics, LinearAlgebra
using Base.Threads
include("../../src/fput_core.jl");   using .FPUTCore
include("../../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../../src/fput_analysis.jl");    using .FPUTAnalysis
include("../../src/block_integration.jl"); using .BlockIntegration

function load_config(path::String)
    d = TOML.parsefile(path)
    phys = d["physics"]
    sim  = d["simulation"]
    out  = d["output"]

    N_values = haskey(phys, "N_values") ? Int.(phys["N_values"]) : [Int(phys["N"])]

    return (
        N_values = N_values,
        boundary = Symbol(phys["boundary"]),
        system_type = phys["system_type"],
        nonlinear = Symbol(phys["nonlinear"]),
        param_values = Float64.(phys["param_values"]),
        delta_values = Float64.(phys["delta_values"]),
        initial_energy = Float64(get(phys, "initial_energy", 0.45)),
        init_mode = Int(get(phys, "init_mode", 1)),
        TMAX = Float64(sim["TMAX"]),
        T_block = Float64(sim["T_block"]),
        DT = Float64(sim["DT"]),
        save_every = Int(sim["save_every"]),
        downsample = Int(sim["downsample"]),
        entropy_delta = Float64(get(sim, "entropy_delta", 0.6)),
        debug = Bool(sim["debug"]),
        base_dir = out["base_dir"],
    )
end

function run_for_N(idx::Int, N::Int, pval::Float64, delta::Float64, cfg)
    println("[N-sweep $idx] Starting N=$N p=$pval Δ=$delta")

    sp = FPUTCore.SystemParams(N, (cfg.system_type=="springs") ? delta : 0.0, (cfg.system_type=="masses") ? delta : 0.0,
                              cfg.nonlinear==:alpha ? pval : 0.0,
                              cfg.nonlinear==:beta ? pval : 0.0,
                              cfg.boundary)

    # Build system and initial condition
    k, m = FPUTCore.make_system(sp)
    q_cur = zeros(N); v_cur = zeros(N)
    freq, V = FPUTCore.find_normal_modes(k, m, sp.boundary)
    target_mode = cfg.init_mode
    amplitude = sqrt(2 * cfg.initial_energy) / freq[target_mode]
    U = Diagonal(1 ./ sqrt.(m)) * V
    q_cur .= amplitude .* U[:, target_mode]

    result = integrate_in_blocks(sp, q_cur, v_cur, N, freq, V, m, target_mode;
                                  TMAX=cfg.TMAX, T_block=cfg.T_block, DT=cfg.DT,
                                  save_every=cfg.save_every, downsample=cfg.downsample,
                                  debug=true, label="N-sweep $idx")

    println("[N-sweep $idx] Done; times=$(length(result.scaled_t)) steps")

    entropy = FPUTAnalysis.spectral_entropy(result.modal_E, cfg.entropy_delta)

    return (N=N, param=pval, Delta=delta, scaled_t=result.scaled_t, modal_E=result.modal_E, entropy=entropy)
end

function main()
    if isempty(ARGS)
        println("Usage: julia scripts/compute_trajectories_Nsweep.jl <config.toml>")
        return
    end
    cfg = load_config(ARGS[1])
    mkpath(cfg.base_dir)

    # Build task list
    tasks = Vector{Tuple{Int,Float64,Float64}}()
    for N in cfg.N_values
        for pval in cfg.param_values
            for delta in cfg.delta_values
                push!(tasks, (N, pval, delta))
            end
        end
    end

    M = length(tasks)
    results = Vector{Any}(undef, M)

    # Parallel execution across tasks (threaded)
    @threads for i in 1:M
        N, pval, delta = tasks[i]
        results[i] = try
            run_for_N(i, N, pval, delta, cfg)
        catch e
            println("[N-sweep $i] ERROR: $e")
            nothing
        end
    end

    valid = filter(!isnothing, results)
    outpath = joinpath(cfg.base_dir, "sweep_N_results_$(Dates.today()).jld2")
    @info "Saving results to $outpath"
    jldsave(outpath; results=valid, config=ARGS[1])
    println("Saved: $outpath with $(length(valid)) entries")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
