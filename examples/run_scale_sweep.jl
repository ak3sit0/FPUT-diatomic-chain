# run_scale_sweep.jl - Versión Simplificada
using TOML, Dates, JLD2

# Configuración básica
TMAX = 1e6
NS = [32, 64, 128, 256]
BASE_CONFIG = "configs/alpha_sweep_N_periodic.toml"
OUTPUT_DIR = "results/raw/alpha_sweep_N_periodic"

# Validar que el config base existe (previene errores silenciosos)
if !isfile(BASE_CONFIG)
    error("Config base no encontrado: $BASE_CONFIG\n" * 
          "Ejecuta este script desde la raíz del proyecto.")
end

println("🚀 Iniciando barrido N=$NS con TMAX=$TMAX")

# 1. Crear configs y lanzar procesos
procs = []
for N in NS
    try
        cfg = TOML.parsefile(BASE_CONFIG)
        cfg["physics"]["N"] = N
        cfg["simulation"]["TMAX"] = TMAX
        
        tmp_path = "configs/tmp_N$(N).toml"
        open(tmp_path, "w") do io; TOML.print(io, cfg); end
        
        # Lanzar sin esperar (wait=false)
        p = run(`julia --project=. examples/compute_modal_energies.jl $tmp_path`, wait=false)
        push!(procs, (p, tmp_path, N))
        println("  ✓ [N=$N] Lanzado correctamente")
    catch e
        println("  ✗ [N=$N] Error: $e")
    end
end

# 2. Esperar y limpiar
if isempty(procs)
    println("⚠️  No hay procesos para ejecutar. Verifica tu configuración.")
    exit(1)
end

println("\n⏳ Esperando a que terminen...")
results_status = []

for (p, path, N) in procs
    try
        wait(p)
        
        if p.exitcode == 0
            # Buscar el archivo .jld2 más reciente en OUTPUT_DIR
            files = filter(f -> endswith(f, ".jld2") && !startswith(basename(f), "tmp_"), 
                          readdir(OUTPUT_DIR; join=true))
            if !isempty(files)
                latest = sort(files; by=mtime)[end]
                # Renombrar para incluir el N: results_springs_alpha_periodic_N32_YYY-MM-DD.jld2
                new_name = replace(basename(latest), ".jld2" => "_N$(N).jld2")
                new_path = joinpath(OUTPUT_DIR, new_name)
                mv(latest, new_path; force=true)
                println("  [N=$N] Terminado: ✓ OK → $(new_name)")
                push!(results_status, true)
            else
                println("  [N=$N] Terminado: ✓ OK (archivo no encontrado)")
                push!(results_status, true)
            end
        else
            println("  [N=$N] Terminado: ✗ Código $(p.exitcode)")
            push!(results_status, false)
        end
    catch e
        println("  [N=$N] Error durante procesamiento: $e")
        push!(results_status, false)
    finally
        rm(path; force=true)
    end
end

n_ok = count(results_status)
println("\n✅ Barrido completado: $n_ok/$(length(procs)) exitosos.")
println("📁 Resultados en: $OUTPUT_DIR/results_springs_alpha_periodic_*.jld2")
