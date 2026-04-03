# run_scale_sweep.jl - Versión Simplificada
using TOML, Dates

# Configuración básica
TMAX = 1e6
NS = [32, 64, 128, 256]
BASE_CONFIG = "configs/alpha_sweep_N_periodic.toml"

println("Iniciando barrido N=$NS con TMAX=$TMAX")

# 1. Crear configs y lanzar procesos
procs = []
for N in NS
    cfg = TOML.parsefile(BASE_CONFIG)
    cfg["physics"]["N"] = N
    cfg["simulation"]["TMAX"] = TMAX
    
    tmp_path = "configs/tmp_N$(N).toml"
    open(tmp_path, "w") do io; TOML.print(io, cfg); end
    
    # Lanzar sin esperar (wait=false)
    p = run(`julia --project=. examples/compute_modal_energies.jl $tmp_path`, wait=false)
    push!(procs, (p, tmp_path, N))
    println("  [N=$N] Lanzado (PID=$(p.pid))")
end

# 2. Esperar y limpiar
println("\n Esperando a que terminen...")
for (p, path, N) in procs
    wait(p)
    rm(path; force=true)
    status = p.exitcode == 0 ? "✓ OK" : "✗ Error"
    println("  [N=$N] Terminado: $status")
end

println("\n Barrido completado.")
