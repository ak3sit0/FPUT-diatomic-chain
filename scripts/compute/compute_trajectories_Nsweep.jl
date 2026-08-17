"""
scripts/compute_trajectories_Nsweep.jl

Barrido sobre N_values definido en un TOML. Guarda un JLD2 con entradas
por cada N: (N, scaled_t, modal_E, entropy).

Usage:
  julia --project=. scripts/compute_trajectories_Nsweep.jl configs/templates/sweep_N_acoustic.toml
"""

using TOML, JLD2, Dates, Statistics, LinearAlgebra
using Base.Threads
include("../src/fput_core.jl");   using .FPUTCore
include("../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../src/fput_analysis.jl");    using .FPUTAnalysis

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

    # accumulators
    T_total = Float64[]
    # collect blocks of modal energies to avoid repeated hcat allocations
    modal_E_blocks = Vector{Matrix{Float64}}()
    t_cur = 0.0

    while t_cur < cfg.TMAX
        t_next = min(t_cur + cfg.T_block, cfg.TMAX)
        saveat = t_cur:cfg.save_every*cfg.DT:t_next

        Qb, Vb, Tb, kvec, mvec = FPUTFastRunner.solve_fput(sp, q_cur, v_cur, (t_cur, t_next), cfg.DT; saveat=saveat)

        # Transform time to cycles using acoustic mode
        Tb = Tb .* freq[target_mode] ./ (2π)

        # Normalize shapes
        if size(Qb,1) == N && size(Qb,2) >= 1
            Qmat = Float64.(Qb)
            Vmat = Float64.(Vb)
        elseif size(Qb,2) == N && size(Qb,1) >= 1
            Qmat = Float64.(Qb')
            Vmat = Float64.(Vb')
        else
            error("Unexpected Q/V shape: Q=$(size(Qb)), V=$(size(Vb))")
        end

        modal_Eb = FPUTAnalysis.compute_modal_energies(Qmat, Vmat, freq, V, m)

        # Append times avoiding duplicate at boundary and applying downsample
        if !isempty(T_total) && !isempty(Tb) && isapprox(T_total[end], Tb[1]; atol=1e-12, rtol=0)
            if length(Tb) > 1
                    idx_ds = 2:cfg.downsample:length(Tb)
                    append!(T_total, Float64.(Tb[idx_ds]))
                    block = Float64.(modal_Eb[:, idx_ds])
                    if size(block, 2) > 0
                        push!(modal_E_blocks, block)
                    end
            end
        else
            idx_ds = 1:cfg.downsample:length(Tb)
            append!(T_total, Float64.(Tb[idx_ds]))
            block = Float64.(modal_Eb[:, idx_ds])
            if size(block, 2) > 0
                push!(modal_E_blocks, block)
            end
        end

        # update state
        if size(Qmat, 2) >= 1
            q_cur .= Qmat[:, end]
            v_cur .= Vmat[:, end]
        end

        t_cur = t_next
        println("[N-sweep $idx] progress: t=$(round(t_cur; digits=3)) / $(cfg.TMAX)")
    end

    println("[N-sweep $idx] Done; times=$(length(T_total)) steps")

    # Combine blocks once to avoid repeated allocations. If no blocks, produce empty matrix.
    if isempty(modal_E_blocks)
        modal_E_total = Array{Float64}(undef, N, 0)
    else
        modal_E_total = reduce(hcat, modal_E_blocks)
    end

    # Compute entropy
    entropy = FPUTAnalysis.spectral_entropy(modal_E_total, cfg.entropy_delta)

    return (N=N, param=pval, Delta=delta, scaled_t=T_total, modal_E=modal_E_total, entropy=entropy)
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
