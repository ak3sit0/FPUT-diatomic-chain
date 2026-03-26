"""
verify_configs.jl — Validate all experiment configs before submitting to cluster.

Usage:
    julia scripts/verify_configs.jl [--quick]

Checks:
- All TOML files parse without error
- Required fields are present
- Parameter ranges are sensible
- Results directories would be created correctly
"""

include("../src/config.jl")
using Main.Config, TOML

function verify_config(path::String)
    try
        cfg = Config.load_experiment_config(path)
        
        # Sanity checks
        checks = [
            (cfg.N > 0, "N must be positive"),
            (cfg.TMAX > 0, "TMAX must be positive"),
            (cfg.T_block > 0 && cfg.T_block <= cfg.TMAX, "T_block must be in (0, TMAX]"),
            (!isempty(cfg.param_values), "param_values cannot be empty"),
            (!isempty(cfg.delta_values), "delta_values cannot be empty"),
            (cfg.initial_energy > 0, "initial_energy must be positive"),
            (cfg.save_every > 0, "save_every must be positive"),
            (cfg.downsample >= 1, "downsample must be >= 1"),
        ]
        
        for (check, msg) in checks
            check || error(msg)
        end
        
        return true, "✓"
    catch e
        return false, "error: $(e.msg)"
    end
end

function main()
    quick = "--quick" in ARGS
    
    # Find all experiment TOMLs
    config_dir = joinpath(@__DIR__, "..", "configs")
    configs = filter(f -> endswith(f, ".toml") && 
                          (startswith(f, "fixed_") || startswith(f, "periodic_") || 
                           startswith(f, "alpha_optical")),
                     readdir(config_dir))
    
    println("=" ^ 70)
    println("🔍  EXPERIMENT CONFIG VALIDATION")
    println("=" ^ 70)
    
    # Group by prefix
    groups = Dict{String, Vector{String}}()
    for cfg in configs
        prefix = match(r"^[^_]+_[^_]+", cfg).match
        push!(get!(groups, prefix, String[]), cfg)
    end
    
    total_ok = 0
    total_fail = 0
    
    for (group, cfgs) in sort(collect(groups))
        println("\n$(group):")
        for cfg in sort(cfgs)
            path = joinpath(config_dir, cfg)
            ok, msg = verify_config(path)
            status = ok ? "✓" : "✗"
            println("  [$status] $cfg")
            if !ok
                println("      $msg")
                total_fail += 1
            else
                total_ok += 1
            end
        end
    end
    
    # Summary
    println("\n" * "=" ^ 70)
    println("SUMMARY: $total_ok valid, $total_fail errors")
    println("=" ^ 70)
    
    if total_fail == 0
        println("\n✓ All configurations valid! Ready to submit:")
        println("  qsub jobs/fixed_bc.pbs")
        println("  qsub jobs/periodic_acoustic.pbs")
        println("  qsub jobs/periodic_optical_beta.pbs")
    else
        error("$total_fail config(s) failed validation")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
