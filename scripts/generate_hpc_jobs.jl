"""
    generate_hpc_jobs.jl

Generador de jobs PBS/Torque para ensambles multi-N en clúster HPC.
Lee una config maestra con `N_values` y genera un .toml + .pbs por cada N.

Ventajas:
  - Un job por N → fallos aislados, relanzables independientemente
  - Recursos escalados automáticamente (cores y walltime según N)
  - Workflow offline: generar → revisar → qsub

Usage:
  julia --project=. scripts/generate_hpc_jobs.jl configs/templates/ensemble_N_sweep.toml

Output:
  - jobs/ensemble_N<N>.toml (config mutada)
  - jobs/ensemble_N<N>.pbs (script PBS)
"""

using TOML
using Dates: today

function resources_for_N(N)
    """Retorna (ppn, walltime) escalados según N."""
    if N <= 64
        return (ppn=8, walltime="12:00:00")
    elseif N <= 128
        return (ppn=16, walltime="24:00:00")
    else  # N = 256
        return (ppn=16, walltime="48:00:00")
    end
end

function write_job_files(base_toml::Dict, N::Int, output_dir::String)
    """Genera .toml y .pbs para un N específico."""
    cfg = deepcopy(base_toml)

    # Remover N_values si existe, poner N fijo
    if haskey(cfg["physics"], "N_values")
        delete!(cfg["physics"], "N_values")
    end
    cfg["physics"]["N"] = N

    # Mutate base_dir con sufijo _N<N>
    original_base = get(cfg["output"], "base_dir", "results/data/ensemble_N_sweep")
    cfg["output"]["base_dir"] = "$(original_base)_N$(N)"

    # Nombres de archivos
    toml_path = joinpath(output_dir, "ensemble_N$(N).toml")
    pbs_path = joinpath(output_dir, "ensemble_N$(N).pbs")

    # Escribir TOML mutado
    open(toml_path, "w") do io
        TOML.print(io, cfg)
    end

    # Recursos escalados
    res = resources_for_N(N)
    ppn = res.ppn
    walltime = res.walltime

    # Escribir PBS script
    open(pbs_path, "w") do io
        println(io, "#!/bin/bash")
        println(io, "#PBS -N ensemble_N$(N)")
        println(io, "#PBS -l nodes=1:ppn=$(ppn)")
        println(io, "#PBS -l walltime=$(walltime)")
        println(io, "#PBS -o results/logs/ensemble_N$(N)_\$(date +%Y%m%d_%H%M%S).log")
        println(io, "#PBS -j oe")
        println(io)
        println(io, "cd \$PBS_O_WORKDIR")
        println(io, "echo \"Starting ensemble N=$(N) on \$(hostname) at \$(date)\"")
        println(io, "echo \"Available cores: \$(nproc)\"")
        println(io)
        println(io, "julia --project=. -t $(ppn) scripts/compute_ensemble.jl $(toml_path)")
        println(io)
        println(io, "echo \"Finished ensemble N=$(N) at \$(date)\"")
    end

    # Resultado
    return (toml=toml_path, pbs=pbs_path, ppn=ppn, walltime=walltime)
end

function main()
    if isempty(ARGS)
        println("Usage: julia scripts/generate_hpc_jobs.jl <config_template.toml>")
        println("Example: julia scripts/generate_hpc_jobs.jl configs/templates/ensemble_N_sweep.toml")
        return
    end

    template_path = ARGS[1]
    if !isfile(template_path)
        error("Config template not found: $template_path")
    end

    base_toml = TOML.parsefile(template_path)
    output_dir = "jobs"
    mkpath(output_dir)
    mkpath("results/logs")

    # Extraer N_values de la config maestra
    N_values = if haskey(base_toml["physics"], "N_values")
        Vector{Int}(base_toml["physics"]["N_values"])
    else
        error("Config debe contener 'N_values' en [physics]")
    end

    println("=== Generador HPC: Ensamble Multi-N ===")
    println("Config: $template_path")
    println("N_values: $N_values")
    println()

    generated = []

    for N in N_values
        result = write_job_files(base_toml, N, output_dir)
        push!(generated, (N=N, result=result))

        println("✓ N=$N")
        println("  TOML:   $(result.toml)")
        println("  PBS:    $(result.pbs)")
        println("  Recurso: ppn=$(result.ppn), walltime=$(result.walltime)")
        println()
    end

    # Sumario y comandos
    println("\n=== Sumario ===")
    println("$(length(generated)) jobs generados en 'jobs/'")
    println()

    println("=== Próximos pasos ===")
    println()
    println("1. Revisar los archivos generados:")
    println("   ls -lh jobs/")
    println()
    println("2. Inspeccionar un PBS (opcional):")
    println("   cat jobs/ensemble_N64.pbs")
    println()
    println("3. Test local (antes de enviar a clúster):")
    println("   julia --project=. -t 4 scripts/compute_ensemble.jl jobs/ensemble_N32.toml")
    println()
    println("4. Enviar todos los jobs al clúster:")
    println("   for job in jobs/ensemble_N*.pbs; do qsub \$job; done")
    println()
    println("5. Monitorear progreso:")
    println("   qstat")
    println("   tail -f results/logs/ensemble_N*.log")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
