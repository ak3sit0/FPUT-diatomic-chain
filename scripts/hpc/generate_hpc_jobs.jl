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
include("../../src/fput_core.jl"); using .FPUTCore

"""
    tmax_for_N(cfg, N) -> TMAX

Si el config trae `scaled_t_max`, deriva TMAX = scaled_t_max·2π/ω₂ para ESTE N.
scaled_t = t·ω₂/2π y ω₂ ∝ 1/N, así que con TMAX FIJO el scaled_t alcanzado cae como
1/N (3.1e4 a N=32 → 3.9e3 a N=256): E_opt/E saldría más bajo a N grande sólo por
haber corrido 8× menos tiempo efectivo, y se leería como "el efecto se debilita".

ω₂ depende de Δκ, y TMAX es único para las 9 Δκ del job ⇒ se toma el TMAX MAYOR
(el de la ω₂ menor) para que TODAS las Δκ alcancen al menos scaled_t_max.
"""
function tmax_for_N(cfg::Dict, N::Int)
    sim = cfg["simulation"]
    haskey(sim, "scaled_t_max") || return Float64(sim["TMAX"])

    phys     = cfg["physics"]
    boundary = Symbol(phys["boundary"])
    stype    = phys["system_type"]
    stmax    = Float64(sim["scaled_t_max"])

    tmaxs = map(Float64.(phys["delta_values"])) do dk
        sp = SystemParams(N, stype == "springs" ? dk : 0.0,
                             stype == "masses"  ? dk : 0.0, 0.0, 0.0, boundary)
        k, m    = make_system(sp)
        freq, _ = find_normal_modes(k, m, boundary)
        sort!(freq)
        stmax * 2π / freq[boundary == :fixed ? 1 : 2]
    end
    maximum(tmaxs)
end

"""
ppn = 20 para que n_real=100 salgan 5 rondas EXACTAS (100/20). Con ppn=16 serían
6.25 → 7 rondas, y la última iría con 4 tareas de 16 (12% de desperdicio).
El paralelismo vive en las realizaciones, así que ppn > n_real no aportaría nada.
"""
function resources_for_N(N)
    if N <= 64
        return (ppn=20, walltime="06:00:00")
    elseif N <= 128
        return (ppn=20, walltime="12:00:00")
    else  # N = 256
        return (ppn=20, walltime="24:00:00")
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

    # save_every = 3N por cada N. Con save_every FIJO, nt = TMAX/(DT·save_every) es
    # constante y modal_E = N·nt crece como N (y como N² si además TMAX ∝ N).
    # Con 3N, nt ∝ 1/N y modal_E queda constante (~13 MB) para cualquier N.
    cfg["simulation"]["save_every"] = 3 * N

    # TMAX ∝ N a scaled_t fijo (ver tmax_for_N). Se resuelve aquí porque el generador
    # es quien conoce N; compute_ensemble.jl sigue leyendo un TMAX absoluto.
    cfg["simulation"]["TMAX"] = tmax_for_N(base_toml, N)
    haskey(cfg["simulation"], "scaled_t_max") && delete!(cfg["simulation"], "scaled_t_max")
    cfg["simulation"]["T_block"] = cfg["simulation"]["TMAX"] / 20

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

    # Escribir PBS script (robusto con error handling)
    log_file = "results/logs/ensemble_N$(N)_\$PBS_JOBID.log"
    open(pbs_path, "w") do io
        println(io, "#!/bin/bash")
        println(io, "#PBS -N ensemble_N$(N)")
        println(io, "#PBS -l nodes=1:ppn=$(ppn)")
        println(io, "#PBS -l walltime=$(walltime)")
        println(io, "#PBS -o $(log_file)")
        println(io, "#PBS -j oe")
        println(io)
        println(io, "cd \$PBS_O_WORKDIR")
        println(io, "mkdir -p results/logs results/data")
        println(io, "LOG=$(log_file)")
        println(io)
        println(io, "{")
        println(io, "  echo \"================================\"")
        println(io, "  echo \"Ensemble N=$(N) - Job ID: \$PBS_JOBID\"")
        println(io, "  echo \"Started at: \$(date '+%Y-%m-%d %H:%M:%S')\"")
        println(io, "  echo \"Hostname: \$(hostname)\"")
        println(io, "  echo \"Available cores: \$(nproc)\"")
        println(io, "  echo \"Running with: $(ppn) threads\"")
        println(io, "  echo \"Walltime limit: $(walltime)\"")
        println(io, "  echo \"================================\"")
        println(io, "  echo")
        println(io)
        println(io, "  # Correr la simulación")
        println(io, "  julia --project=. -t $(ppn) scripts/compute_ensemble.jl $(toml_path)")
        println(io, "  EXIT_CODE=\$?")
        println(io)
        println(io, "  # Verificar que se guardó el archivo de resultados")
        println(io, "  RESULT_FILE=\$(find $(cfg["output"]["base_dir"]) -name '*.jld2' 2>/dev/null | head -1)")
        println(io, "  if [ -z \"\$RESULT_FILE\" ] && [ \$EXIT_CODE -eq 0 ]; then")
        println(io, "    echo \"⚠ WARNING: Script completó pero NO se encontró archivo .jld2\"")
        println(io, "    ls -lh $(cfg["output"]["base_dir"])/ || echo \"Directorio no existe\"")
        println(io, "    EXIT_CODE=1")
        println(io, "  elif [ ! -z \"\$RESULT_FILE\" ]; then")
        println(io, "    echo \"✓ Archivo de resultados: \$RESULT_FILE\"")
        println(io, "    ls -lh \$RESULT_FILE")
        println(io, "  fi")
        println(io)
        println(io, "  echo")
        println(io, "  echo \"================================\"")
        println(io, "  if [ \$EXIT_CODE -eq 0 ]; then")
        println(io, "    echo \"✓ JOB COMPLETED SUCCESSFULLY\"")
        println(io, "    echo \"Finished at: \$(date '+%Y-%m-%d %H:%M:%S')\"")
        println(io, "    touch results/logs/ensemble_N$(N)_\$PBS_JOBID.SUCCESS")
        println(io, "  else")
        println(io, "    echo \"✗ JOB FAILED with exit code \$EXIT_CODE\"")
        println(io, "    echo \"Failed at: \$(date '+%Y-%m-%d %H:%M:%S')\"")
        println(io, "    touch results/logs/ensemble_N$(N)_\$PBS_JOBID.FAILED")
        println(io, "  fi")
        println(io, "  echo \"================================\"")
        println(io, "} 2>&1 | tee -a \$LOG")
        println(io)
        println(io, "exit \$EXIT_CODE")
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
