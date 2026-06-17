"""
    fput_simulation.jl

A Literate Programming rewrite of the FPUT simulation engine.

Layers:
1. Configuration & Data Structures
2. Physics Core (FPUTCore)
3. Integration & Running (FPUTFastRunner)
4. Analysis & Metrics (FPUTAnalysis)
5. Orchestration (Main)

Usage:
    julia --project=. examples/fput_simulation.jl [config.toml]
"""

# ── Dependencies ──
using TOML, JLD2, Dates, LinearAlgebra, Random
include("../src/fput_core.jl");        using .FPUTCore
include("../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../src/fput_analysis.jl");    using .FPUTAnalysis

# ============================================================================
# Layer 1: Configuration Loading
# ============================================================================

function load_clean_config(path::String)
    d = TOML.parsefile(path)
    p = d["physics"]
    s = d["simulation"]
    
    # Map TOML to SystemParams
    params = SystemParams(
        Int(p["N"]),
        Float64(get(p, "DeltaK", 0.0)),
        Float64(get(p, "DeltaM", 0.0)),
        (Symbol(p["nonlinear"]) == :alpha) ? Float64(p["param_values"][1]) : 0.0,
        (Symbol(p["nonlinear"]) == :beta)  ? Float64(p["param_values"][1]) : 0.0,
        Symbol(p["boundary"])
    )
    
    return params, d["physics"], s, d["output"]["base_dir"]
end

# ============================================================================
# Layer 2: Initial Condition Strategy
# ============================================================================

function prepare_initial_condition(p::SystemParams, k, m, phys_meta, sim_meta)
    freq, V = find_normal_modes(k, m, p.boundary)
    U = Diagonal(1 ./ sqrt.(m)) * V # Normalization matrix
    
    energy = Float64(get(phys_meta, "initial_energy", 0.45))
    init_type = get(phys_meta, "init_type", "mode")
    
    q0 = zeros(p.N)
    v0 = zeros(p.N)
    
    if init_type == "mode"
        target_mode = haskey(phys_meta, "init_mode") ? Int(phys_meta["init_mode"]) : (p.boundary == :fixed ? 1 : 2)
        println("  -> Initializing specific mode: $target_mode with Energy = $energy")
        # E = 0.5 * omega^2 * A^2
        amplitude = sqrt(2 * energy) / freq[target_mode]
        q0 .= amplitude .* U[:, target_mode]
        
    elseif init_type == "random"
        println("  -> Initializing random thermal state with Total Energy ≈ $energy")
        # Distribute energy equally among modes (equipartition starting attempt)
        # Random phases in [0, 2π)
        for i in 1:p.N
            omega = freq[i]
            if omega > 1e-10 # skip zero modes
                A = sqrt(2 * (energy / p.N)) / omega
                phi = rand() * 2π
                # q_modal = A * cos(phi), v_modal = -A * omega * sin(phi)
                q_modal = A * cos(phi)
                v_modal = -A * omega * sin(phi)
                
                q0 .+= q_modal .* U[:, i]
                v0 .+= v_modal .* U[:, i]
            end
        end
    else
        error("Unknown init_type: $init_type. Expected 'mode' or 'random'.")
    end
    
    return q0, v0, freq, V
end

# ============================================================================
# Layer 5: Main Entry Point (The Narrative)
# ============================================================================

function run_fput_experiment(config_path::String)
    println("--- Starting FPUT Refactored Simulation ---")
    
    # 1. Load context
    sys_params, phys_meta, sim_meta, base_dir = load_clean_config(config_path)
    k, m = make_system(sys_params)
    
    # 2. Setup Initial Condition
    q0, v0, freqs, V_modes = prepare_initial_condition(sys_params, k, m, phys_meta, sim_meta)

    
    # 3. Solve the Trajectory
    tspan = (0.0, Float64(sim_meta["TMAX"]))
    dt = get(sim_meta, "DT", 0.05)
    saveat = tspan[1]:(sim_meta["save_every"] * dt):tspan[2]
    
    println("Solving ODE Grid (N=$(sys_params.N), TMAX=$(tspan[2]))...")
    Q, V, T, _, _ = solve_fput(sys_params, q0, v0, tspan, dt; saveat=saveat)
    
    # 4. Perform Analysis Pipeline
    println("Analyzing Energies and Entropy...")
    modal_E = compute_modal_energies(Q, V, freqs, V_modes, m)
    entropy = spectral_entropy(modal_E, 0.6)
    
    # 5. Guardar resultados
    outpath = joinpath(base_dir, "refactored_results_$(Dates.today()).jld2")
    mkpath(base_dir)
    # Save primitive config dict to avoid type reconstruction warnings on load
    cfg_map = Dict(
        "N" => sys_params.N,
        "DeltaK" => sys_params.delta_k,
        "DeltaM" => sys_params.delta_m,
        "alpha" => sys_params.alpha,
        "beta" => sys_params.beta,
        "boundary" => String(sys_params.boundary),
    )
    jldsave(outpath; modal_E=modal_E, entropy=entropy, time=T, config=cfg_map)
    
    println("Completed. Result saved in: $outpath")
    println("Final Entropy: $(round(entropy[end], digits=4))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    cfg_file = isempty(ARGS) ? "configs/templates/test_quick.toml" : ARGS[1]
    run_fput_experiment(cfg_file)
end
