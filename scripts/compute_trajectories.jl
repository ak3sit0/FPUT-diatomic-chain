"""
    compute_trajectories.jl

Simple, robust reimplementation of `compute_modal_energies.jl` using the
refactored `src/` modules. Performs parameter sweeps (param_values × delta_values),
saves a single `.jld2` with multiple run entries and minimal reproducibility metadata.

Usage:
  julia --project=. scripts/compute_trajectories.jl configs/templates/test_quick.toml
"""

using TOML, JLD2, Dates, Base.Threads, Statistics, LinearAlgebra
include("../src/fput_core.jl");   using .FPUTCore
include("../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../src/fput_analysis.jl");    using .FPUTAnalysis

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

    # Prepare accumulators
    T_total = Float64[]
    # collect blocks of modal energies to avoid repeated hcat allocations
    modal_E_blocks = Vector{Matrix{Float64}}()
    t_cur = 0.0

    # Integrate in blocks so we can print progress like the recovery script
    while t_cur < cfg.TMAX
        t_next = min(t_cur + cfg.T_block, cfg.TMAX)
        saveat = t_cur:cfg.save_every*cfg.DT:t_next

        Qb, Vb, Tb, kvec, mvec = FPUTFastRunner.solve_fput(sp, q_cur, v_cur, (t_cur, t_next), cfg.DT; saveat=saveat)

        # Transform time to cycles of oscillation (normalized by acoustic mode frequency)
        Tb = Tb .* freq[target_mode] ./ (2π)

        # Normalize shapes: ensure Qmat, Vmat are (N × nt)
        if size(Qb,1) == cfg.N && size(Qb,2) >= 1
            Qmat = Float64.(Qb)
            Vmat = Float64.(Vb)
        elseif size(Qb,2) == cfg.N && size(Qb,1) >= 1
            Qmat = Float64.(Qb')
            Vmat = Float64.(Vb')
        else
            error("Unexpected Q/V shape: Q=$(size(Qb)), V=$(size(Vb))")
        end

        # Compute modal energies for this block (modal_Eb is N × nt)
        modal_Eb = FPUTAnalysis.compute_modal_energies(Qmat, Vmat, freq, V, m)

        # Append times avoiding duplicate at boundary
        if !isempty(T_total) && !isempty(Tb) && isapprox(T_total[end], Tb[1]; atol=1e-12, rtol=0)
            # drop the duplicated first column and apply downsample
                if length(Tb) > 1
                idx_ds = 2:cfg.downsample:length(Tb)
                append!(T_total, Float64.(Tb[idx_ds]))
                block = Float64.(modal_Eb[:, idx_ds])
                if size(block,2) > 0
                    push!(modal_E_blocks, block)
                end
            end
        else
            # initial or non-overlapping: apply downsample
            idx_ds = 1:cfg.downsample:length(Tb)
            append!(T_total, Float64.(Tb[idx_ds]))
            block = Float64.(modal_Eb[:, idx_ds])
            if size(block,2) > 0
                push!(modal_E_blocks, block)
            end
        end

        # Update current state to last sample for next block
        if size(Qmat, 2) >= 1
            q_cur .= Qmat[:, end]
            v_cur .= Vmat[:, end]
        end

        t_cur = t_next
        println("[case $idx] progress: t=$(round(t_cur; digits=3)) / $(cfg.TMAX)")
    end

    

    println("[case $idx] Done; times=$(length(T_total)) steps")

    if isempty(modal_E_blocks)
        modal_E_total = Array{Float64}(undef, cfg.N, 0)
    else
        modal_E_total = reduce(hcat, modal_E_blocks)
    end

    return (param=pval, Delta=delta, scaled_t=T_total, modal_E=modal_E_total)
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
