# Compute modal energies: Functional, TOML-driven version
# Usage: julia --project=. --threads=4 examples/compute_modal_energies.jl configs/my_experiment.toml

include("../src/config.jl");           using Main.Config
include("../src/parameters.jl");       using Main.Parameters
include("../src/simulation_runner.jl"); using Main.SimulationRunner
include("../src/energy_analysis.jl");   using Main.EnergyAnalysis

using Base.Threads, LinearAlgebra, JLD2, Dates, Statistics

# ============================================================================
# UTILITIES & HELPERS
# ============================================================================

# Context manager for logging
function with_logger(f::Function, dir::String)
    path = joinpath(dir, "log_$(Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")).txt")
    open(path, "w") do io
        log_fn(msg) = (println(msg); write(io, "$msg\n"); flush(io))
        f(log_fn)
    end
end

# Compact Validation & Stats
function analyze_block_health(modal_E, last_E_tot)
    flat = vec(modal_E)
    valid_flat = filter(!isnan, flat)
    
    stats = (total=sum(valid_flat), mean=mean(valid_flat), std=std(valid_flat),
             n_nan=count(isnan, flat), n_inf=count(isinf, flat))
    
    warnings = String[]
    stats.n_nan > 0 && push!(warnings, "NaNs detected: $(stats.n_nan)")
    stats.n_inf > 0 && push!(warnings, "Infs detected: $(stats.n_inf)")
    stats.total > 1e10 && push!(warnings, "High energy: $(stats.total)")
    
    # Continuity check
    start_E, end_E = sum(modal_E[:, 1]), sum(modal_E[:, end])
    if !isnothing(last_E_tot) && abs(start_E - last_E_tot) / (last_E_tot + 1e-10) > 0.01
        push!(warnings, "Energy jump > 1%")
    end
    
    return stats, warnings, end_E
end

# Compression logic
function compress(modal_E, t, factor)
    idx = 1:factor:size(modal_E, 2)
    (Float32.(modal_E[:, idx]), Float32.(t[idx]), sizeof(modal_E[:, idx])/1e6)
end

# ============================================================================
# CONFIGURATION & STATE MANAGEMENT
# ============================================================================

function build_config(config_path::String)
    cfg    = Config.load_experiment_config(config_path)
    outdir = Config.config_outdir(cfg)
    mkpath(joinpath(outdir, "checkpoints"))
    mkpath(joinpath(outdir, "debug"))

    return (
        system_type    = cfg.system_type,
        nonlinear      = cfg.nonlinear,
        initial        = cfg.initial_condition,   # :low | :high
        boundary       = cfg.boundary,
        N              = cfg.N,
        TMAX           = cfg.TMAX,
        T_block        = cfg.T_block,
        save_every     = cfg.save_every,
        downsample     = cfg.downsample,
        debug          = cfg.debug,
        chk_dir        = joinpath(outdir, "checkpoints"),
        debug_dir      = joinpath(outdir, "debug"),
        outdir         = outdir,
        params         = cfg.param_values,
        deltas         = cfg.delta_values,
        initial_energy = cfg.initial_energy,
        init_mode      = cfg.init_mode,   # nothing → derive from initial_condition
    )
end

function load_or_init_state(chkfile, N::Int)
    if isfile(chkfile)
        return load(chkfile)
    else
        return Dict{String, Any}(
            "q0"        => nothing,
            "v0"        => nothing,
            "t_cur"     => 0.0,
            "last_E"    => nothing,
            "scaled_ts" => Float32[],
            "modal_E"   => Array{Float32}(undef, N, 0),
        )
    end
end

function save_checkpoint(chkfile, state, param, Delta, config)
    jldsave(chkfile; 
        Dict(Symbol(k)=>v for (k,v) in state)..., 
        param=param, 
        Delta=Delta, 
        nonlinear=config.nonlinear, 
        system_type=config.system_type
    )
end

# ============================================================================
# CORE PHYSICS ENGINE
# ============================================================================

function simulate_single_block(state, config, param, Delta)
    t_next = min(state["t_cur"] + config.T_block, config.TMAX)
    
    # Parameter mapping
    alpha = config.nonlinear == :alpha ? param : 0.0
    beta  = config.nonlinear == :beta  ? param : 0.0
    DeltaK = config.system_type == :springs ? Delta : 0.0
    DeltaM = config.system_type == :masses ? Delta : 0.0
    
    # Initial mode: explicit TOML value > initial_condition keyword > defaults
    init_mode = if !isnothing(config.init_mode)
        config.init_mode                             # explicit numeric mode from TOML
    elseif config.initial == :high
        config.N ÷ 2 + 1                             # lowest optical mode (first above gap)
    else
        config.boundary == :fixed ? 1 : 2            # lowest acoustic mode
    end

    # Run one block of integration
    modal_E_block, t_blk, _, v_fin, q_fin = SimulationRunner.run_and_analyze(
        N=config.N, TMAX=t_next, DT=Parameters.DT,
        alpha=alpha, beta=beta, DeltaK=DeltaK, DeltaM=DeltaM,
        boundary=config.boundary, init_mode=init_mode,
        q0_override    = state["q0"],
        v0_override    = state["v0"],
        energy         = config.initial_energy,
        tspan_override = (state["t_cur"], t_next),
        save_every     = config.save_every,
    )

    # Analyze & Compress
    stats, warns, last_E = analyze_block_health(modal_E_block, state["last_E"])
    E_c, t_c, size_mb = compress(modal_E_block, t_blk, config.downsample)
    
    # Update State
    state["modal_E"] = hcat(state["modal_E"], E_c)
    append!(state["scaled_ts"], t_c)
    state["q0"] = q_fin
    state["v0"] = v_fin
    state["t_cur"] = t_next
    state["last_E"] = last_E
    
    return state, stats, warns, size_mb, modal_E_block
end

function simulate_trajectory(idx, param, Delta, config, logger)
    chkfile = joinpath(config.chk_dir, "chk_$(idx).jld2")
    
    logger("[$(threadid())] ↻ Initializing/Resuming $idx (p=$param, d=$Delta)")
    state = load_or_init_state(chkfile, config.N)
    
    try
        while state["t_cur"] < config.TMAX
            # 1. Advance physics by one block
            state, stats, warns, size_mb, modal_E_block = simulate_single_block(state, config, param, Delta)
            
            # 2. Save progress
            save_checkpoint(chkfile, state, param, Delta, config)
            
            # 3. Report status
            status = isempty(warns) ? "✓" : " $(join(warns, "; "))"
            logger("[$(threadid())] $status i=$idx t=$(round(state["t_cur"])) | Mem: $(round(size_mb, digits=1))MB | E_mean: $(round(stats.mean, sigdigits=3))")

            # 4. Handle bad blocks
            if !isempty(warns) && config.debug
                @save joinpath(config.debug_dir, "bad_block_$(idx)_t$(round(Int, state["t_cur"])).jld2") modal_E=modal_E_block stats warns
            end
            
            # 5. Graceful Exit
            isfile(joinpath(config.outdir, "STOP")) && break
        end
        
        # Cleanup & Return
        rm(chkfile; force=true)
        logger("[$(threadid())]  DONE: p=$param, d=$Delta")
        return (param=param, Delta=Delta, scaled_t=state["scaled_ts"], modal_E=state["modal_E"])

    catch e
        logger("[$(threadid())]  ERROR i=$idx: $e")
        config.debug && @save joinpath(config.debug_dir, "crash_$(idx).jld2") e param Delta
        return nothing
    end
end

# ============================================================================
# MAIN ENTRY
# ============================================================================

function main()
    if isempty(ARGS)
        error("Usage: julia compute_modal_energies.jl <config.toml>\n" *
              "Example: julia --project=. --threads=4 " *
              "examples/compute_modal_energies.jl configs/alpha_sweep_periodic.toml")
    end

    config_path = ARGS[1]
    config = build_config(config_path)
    pairs  = vec([(p, d) for p in config.params, d in config.deltas])

    with_logger(config.outdir) do logger
        logger("STARTING: $(length(pairs)) sims | TMAX=$(config.TMAX) | " *
               "Block=$(config.T_block) | Config: $config_path")

        # Parallel execution
        results = Vector{Any}(undef, length(pairs))
        Threads.@threads for i in eachindex(pairs)
            results[i] = simulate_trajectory(
                i, pairs[i][1], pairs[i][2], config, logger)
        end

        # Save results + full reproducibility metadata
        final_file = joinpath(config.outdir,
            "results_$(config.system_type)_$(config.nonlinear)_" *
            "$(config.boundary)_$(Dates.today()).jld2")
        valid_res = filter(!isnothing, results)
        Config.save_with_metadata(final_file, config_path; results=valid_res, config=config)
        logger("COMPLETED. Saved $(length(valid_res))/$(length(pairs)) to $final_file")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
