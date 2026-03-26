# Generate a PBS job file from an experiment TOML config.
# Usage: julia scripts/generate_pbs.jl configs/my_experiment.toml [options]
#
# Options:
#   --ppn=N           Processors per node    (default: 8)
#   --nodes=N         Number of nodes        (default: 1)
#   --mem=Xgb         Memory                 (default: 128gb)
#   --walltime=H:M:S  Wall time              (default: 24:00:00)
#   --queue=NAME      Queue name             (default: workq)

using TOML

function generate_pbs(config_path::String;
                      ppn::Int       = 8,
                      nodes::Int     = 1,
                      mem::String    = "128gb",
                      walltime::String = "24:00:00",
                      queue::String  = "workq")

    isfile(config_path) || error("Config file not found: $config_path")
    d    = TOML.parsefile(config_path)
    name = d["experiment"]["name"]

    # Escape the config path for inline use inside the heredoc
    abs_config = abspath(config_path)

    pbs = """
#!/bin/bash
#PBS -N $(name)
#PBS -q $(queue)
#PBS -l nodes=$(nodes):ppn=$(ppn)
#PBS -l mem=$(mem)
#PBS -l walltime=$(walltime)
#PBS -o results/logs/\${PBS_JOBID}_$(name).log
#PBS -j oe

set -euo pipefail
cd "\$PBS_O_WORKDIR"

JOBTAG="\${PBS_JOBID:-\$(date +%s)}"
export JULIA_NUM_THREADS=$(ppn)
JULIA_BIN="\${JULIA_BIN:-\$HOME/.juliaup/bin/julia}"

echo "Job [\$JOBTAG] started on \$(hostname) at \$(date)"
echo "Config: $(config_path)"

mkdir -p results/logs

"\$JULIA_BIN" --project=. --threads=\$JULIA_NUM_THREADS \\
    examples/compute_modal_energies.jl $(config_path)

RC=\$?
echo "Job finished with code: \$RC at \$(date)"
exit \$RC
"""

    mkpath("jobs")
    outpath = joinpath("jobs", "$(name).pbs")
    write(outpath, lstrip(pbs))
    println("Generated: $outpath")
    return outpath
end

function main()
    if isempty(ARGS)
        error("Usage: julia scripts/generate_pbs.jl configs/my_experiment.toml " *
              "[--ppn=N] [--nodes=N] [--mem=Xgb] [--walltime=H:M:S] [--queue=NAME]")
    end

    config_path = ARGS[1]

    # Parse optional --key=value flags
    kwargs = Dict{Symbol, Any}()
    for arg in ARGS[2:end]
        m = match(r"^--(\w+)=(.+)$", arg)
        isnothing(m) && continue
        key, raw = Symbol(m[1]), m[2]
        kwargs[key] = key in (:ppn, :nodes) ? parse(Int, raw) : String(raw)
    end

    generate_pbs(config_path; kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
