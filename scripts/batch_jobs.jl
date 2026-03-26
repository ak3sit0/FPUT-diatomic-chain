"""
batch_jobs.jl — Generate a single PBS job that runs multiple simulation configs sequentially.

Usage:
    julia scripts/batch_jobs.jl jobs/group_fixed_bc.toml [--ppn=16] [--walltime=HH:MM:SS] ...

Reads a job group TOML file listing multiple experiment configs, then generates a single PBS
that sequentially runs all configs using the parent julia instance (not forked).
"""

using TOML

function generate_batch_pbs(group_path::String;
                           ppn::Int=16,
                           nodes::Int=1,
                           mem::String="128gb",
                           walltime::String="72:00:00",
                           queue::String="workq")
    
    isfile(group_path) || error("Group file not found: $group_path")
    
    group = TOML.parsefile(group_path)
    job_name = group["job"]["name"]
    sims = group["simulations"]
    
    # Validate all configs exist
    for sim in sims
        config_path = sim["config"]
        isfile(config_path) || error("Config not found: $config_path")
    end
    
    outpath = joinpath("jobs", "$(job_name).pbs")
    
    # Build list of commands (one per config, chained with &&)
    config_list = [sim["config"] for sim in sims]
    run_cmds = join([
        "\"\$JULIA_BIN\" --project=. --threads=\$JULIA_NUM_THREADS examples/compute_modal_energies.jl $cfg"
        for cfg in config_list
    ], " && \\\n    ")
    
    pbs_content = """#!/bin/bash
#PBS -N $(job_name)
#PBS -q $queue
#PBS -l nodes=$nodes:ppn=$ppn
#PBS -l mem=$mem
#PBS -l walltime=$walltime
#PBS -o results/logs/\${PBS_JOBID}_$(job_name).log
#PBS -j oe

set -euo pipefail
cd "\$PBS_O_WORKDIR"

JOBTAG="\${PBS_JOBID:-\$(date +%s)}"
export JULIA_NUM_THREADS=$ppn
JULIA_BIN="\${JULIA_BIN:-\$HOME/.juliaup/bin/julia}"

echo "Job [\$JOBTAG] started on \$(hostname) at \$(date)"
echo "Job group: $group_path"
echo "Configs to run: $(length(config_list))"
for cfg in $(join(config_list, " ")); do
    echo "  - \$cfg"
done
echo ""

mkdir -p results/logs

# Run all configs sequentially
$run_cmds

RC=\$?
echo "Job finished with code: \$RC at \$(date)"
exit \$RC
"""
    
    open(outpath, "w") do f
        write(f, pbs_content)
    end
    
    println("Generated: $outpath")
    println("  Configs: $(join(config_list, ", "))")
    println("  PPNs: $ppn, Walltime: $walltime")
    return outpath
end

function main()
    if isempty(ARGS)
        error("Usage: julia scripts/batch_jobs.jl jobs/group_NAME.toml " *
              "[--ppn=N] [--nodes=N] [--mem=Xgb] [--walltime=H:M:S] [--queue=NAME]")
    end

    group_path = ARGS[1]

    # Parse optional --key=value flags
    kwargs = Dict{Symbol, Any}()
    for arg in ARGS[2:end]
        m = match(r"^--(\w+)=(.+)$", arg)
        isnothing(m) && continue
        key, raw = Symbol(m[1]), m[2]
        kwargs[key] = key in (:ppn, :nodes) ? parse(Int, raw) : String(raw)
    end

    generate_batch_pbs(group_path; kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
