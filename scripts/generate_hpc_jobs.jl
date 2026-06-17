"""
    generate_hpc_jobs.jl

Orquestador para barridos de parámetros en clúster (HPC).
Genera archivos .pbs y .toml para cada configuración sin ejecutarlos, 
permitiendo al usuario encolar (qsub) de manera offline e independiente.

Usage:
    julia scripts/generate_hpc_jobs.jl --template configs/templates/test_quick.toml --sweep N
"""

using TOML

function write_job(base_toml, sweep_name, param_name, param_val)
    # Genera una copia mutada del TOML
    cfg = deepcopy(base_toml)
    
    # Mutar el parámetro según la semántica
    if param_name == "N"
        cfg["physics"]["N"] = param_val
    elseif param_name == "alpha"
        cfg["physics"]["nonlinear"] = "alpha"
        cfg["physics"]["param_values"] = [param_val]
    elseif param_name == "beta"
        cfg["physics"]["nonlinear"] = "beta"
        cfg["physics"]["param_values"] = [param_val]
    end
    
    # Nombramiento
    job_id = "$(sweep_name)_$(param_name)_$(param_val)"
    toml_path = "jobs/$job_id.toml"
    pbs_path = "jobs/$job_id.pbs"
    
    # Escribir el nuevo TOML
    open(toml_path, "w") do io
        TOML.print(io, cfg)
    end
    
    # Escribir la macro PBS línea por línea (evita problemas de parsing con interpolación)
    open(pbs_path, "w") do io
        println(io, "#!/bin/bash")
        println(io, "#PBS -N FPUT_$(job_id)")
        println(io, "#PBS -l nodes=1:ppn=8")
        println(io, "#PBS -l walltime=24:00:00")
        println(io, "#PBS -o results/raw/log_$(job_id).out")
        println(io, "#PBS -j oe")
        println(io)
        println(io, "cd \$PBS_O_WORKDIR")
        println(io, "julia --project=. examples/fput_literate_sim.jl $toml_path")
    end
    
    return pbs_path
end

function generate_sweep(template_path::String, sweep_param::String)
    println("--- Generando Sweep HPC para $sweep_param ---")
    mkpath("jobs")
    
    base_toml = TOML.parsefile(template_path)
    
    # Array de valores hipotéticos para el barrido
    sweep_values = if sweep_param == "N"
        [32, 64, 128, 256, 512]
    elseif sweep_param == "alpha"
        [0.05, 0.1, 0.25, 0.5, 1.0]
    elseif sweep_param == "beta"
        [0.05, 0.1, 0.25, 0.5, 1.0]
    else
        error("Sweep parametro no soportado: $sweep_param (Use N, alpha o beta)")
    end
    
    generated = []
    
    for val in sweep_values
        pbs = write_job(base_toml, "HPC_Sweep", sweep_param, val)
        push!(generated, pbs)
        println("  · Generado: $pbs")
    end
    
    println("\\nPara lanzar los trabajos, usa:")
    println("  for job in jobs/*.pbs; do qsub \\\$job; done")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) >= 2 && ARGS[1] == "--template"
        template = ARGS[2]
        sweep_param = length(ARGS) >= 4 && ARGS[3] == "--sweep" ? ARGS[4] : "N"
        generate_sweep(template, sweep_param)
    else
        println("Uso: julia scripts/generate_hpc_jobs.jl --template <ruta.toml> --sweep <N|alpha|beta>")
    end
end
