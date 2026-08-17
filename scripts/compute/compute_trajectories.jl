"""
    compute_trajectories.jl

Simple, robust reimplementation of `compute_modal_energies.jl` using the
refactored `src/` modules. Performs parameter sweeps (param_values × delta_values),
saves a single `.jld2` with multiple run entries and minimal reproducibility metadata.

Usage:
  julia --project=. scripts/compute_trajectories.jl configs/templates/test_quick.toml
"""

using TOML, JLD2, Dates, Base.Threads, Statistics, LinearAlgebra
include("../../src/fput_core.jl");   using .FPUTCore
include("../../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../../src/fput_analysis.jl");    using .FPUTAnalysis
include("../../src/block_integration.jl"); using .BlockIntegration

function build_config(path::String)
    d = TOML.parsefile(path)
    phys = d["physics"]
    sim  = d["simulation"]
    out  = d["output"]

    return (
        N = Int(phys["N"]),
        boundary = Symbol(phys["boundary"]),
        system_type = phys["system_type"],
        nonlinear = Symbol(phys["nonlinear"]),
        param_values = Float64.(phys["param_values"]),
        delta_values = Float64.(phys["delta_values"]),
        initial_energy = Float64(get(phys, "initial_energy", 0.45)),
        init_type = get(phys, "initial_condition", "low"),
        TMAX = Float64(sim["TMAX"]),
        T_block = Float64(sim["T_block"]),
        DT = Float64(sim["DT"]),
        save_every = Int(sim["save_every"]),
        downsample = Int(sim["downsample"]),
        debug = Bool(sim["debug"]),
        base_dir = out["base_dir"],
    )
end

function run_case(idx, pval, delta, cfg)
    println("[case $idx] Starting p=$pval Δ=$delta")
    # Build SystemParams
    sp = FPUTCore.SystemParams(cfg.N, (cfg.system_type=="springs") ? delta : 0.0, (cfg.system_type=="masses") ? delta : 0.0,
                              cfg.nonlinear==:alpha ? pval : 0.0,
                              cfg.nonlinear==:beta ? pval : 0.0,
                              cfg.boundary)

    # Initial condition: mode or random (simple: mode)
    k, m = FPUTCore.make_system(sp)
    q_cur = zeros(cfg.N); v_cur = zeros(cfg.N)
    freq, V = FPUTCore.find_normal_modes(k, m, sp.boundary)

    # Sorting modes by frequency to ensure consistent mode selection across cases
    idx_sort = sortperm(freq)
    freq = freq[idx_sort]
    V    = V[:, idx_sort]

    # default to lowest acoustic mode
    target_mode = sp.boundary == :fixed ? 1 : 2
    amplitude = sqrt(2 * cfg.initial_energy) / freq[target_mode]
    U = Diagonal(1 ./ sqrt.(m)) * V
    q_cur .= amplitude .* U[:, target_mode]

    result = integrate_in_blocks(sp, q_cur, v_cur, cfg.N, freq, V, m, target_mode;
                                  TMAX=cfg.TMAX, T_block=cfg.T_block, DT=cfg.DT,
                                  save_every=cfg.save_every, downsample=cfg.downsample,
                                  debug=true, label="case $idx")

    println("[case $idx] Done; times=$(length(result.scaled_t)) steps")

    return (param=pval, Delta=delta, scaled_t=result.scaled_t, modal_E=result.modal_E)
end

function main()
    if isempty(ARGS)
        println("Usage: julia scripts/compute_trajectories.jl <config.toml>")
        return
    end
    cfg = build_config(ARGS[1])
    mkpath(cfg.base_dir)

    pairs = collect(Iterators.product(cfg.param_values, cfg.delta_values))
    results = Vector{Any}(undef, length(pairs))

    Threads.@threads for i in eachindex(pairs)
        pval, delta = pairs[i]
        results[i] = try
            run_case(i, pval, delta, cfg)
        catch e
            println("[case $i] ERROR: $e")
            nothing
        end
    end

    valid = filter(!isnothing, results)
    outpath = joinpath(cfg.base_dir, "sweep_results_$(Dates.today()).jld2")
    @info "Saving results to $outpath"
    jldsave(outpath; results=valid, config=ARGS[1])
    println("Saved: $outpath with $(length(valid)) entries")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
