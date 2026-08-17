"""
    generate_pbs_nsweep.jl

Genera un .toml + .pbs por cada N de un barrido nsweep, con el walltime DERIVADO
del modelo de coste medido en lugar de estimado a ojo.

Uso:
  julia --project=. scripts/generate_pbs_nsweep.jl configs/cases/nsweep_rise_N256.toml
  julia --project=. scripts/generate_pbs_nsweep.jl <config> --ns-per-step-site=95 --safety=2.5

Opciones:
  --ns-per-step-site=C   Constante de la máquina (por defecto 67, medida en el portátil
                         tras la reescritura del kernel). MEDIRLA EN EL CLÚSTER con
                         scripts/calibrate_cost.jl — los nodos suelen ser más lentos.
  --safety=F             Margen sobre el walltime estimado (por defecto 2.5).
  --outdir=DIR           Dónde escribir los jobs (por defecto "jobs").

Modelo de coste (todo verificado por medición):
  ω₂ se calcula exactamente por (N, Δκ)  —NO se aproxima: el ω₂·N ≈ 6.25 sólo vale
  para Δκ pequeño, y la rama acústica se ablanda mucho al crecer Δκ.
  TMAX      = scaled_t_max · 2π/ω₂          (⇒ TMAX ∝ N a scaled_t fijo)
  t_caso    = C · 1e-9 · N · TMAX/dt        (⇒ coste ∝ N² a scaled_t fijo)
  walltime  = t_caso · rondas · contención · safety
  con rondas = ceil(nΔκ / ppn) y contención ≈ 1.7 medida con 4 casos concurrentes.

Un job por N (no por (N,Δκ)): los Δκ del mismo N cuestan lo mismo, así que llenan los
hilos de forma pareja. Una sola trayectoria NO se puede paralelizar en el tiempo, así
que ppn > nΔκ no aporta nada.
"""

using TOML, Printf
include("../../src/fput_core.jl"); using .FPUTCore

const CONTENTION = 1.7   # medido: 4 casos concurrentes vs 1 aislado

getopt(args, key, default) = begin
    i = findfirst(a -> startswith(a, "--$key="), args)
    isnothing(i) ? default : split(args[i], "=", limit=2)[2]
end

"ω del modo de referencia (modo 2 en PBC, modo 1 en frontera fija), exacto."
function omega_ref(N::Int, delta::Float64, boundary::Symbol, system_type::String)
    sp = SystemParams(N,
                      system_type == "springs" ? delta : 0.0,
                      system_type == "masses"  ? delta : 0.0,
                      0.0, 0.0, boundary)
    k, m = make_system(sp)
    freq, _ = find_normal_modes(k, m, boundary)
    sort!(freq)
    freq[boundary == :fixed ? 1 : 2]
end

hhmmss(s) = (h = floor(Int, s/3600); m = floor(Int, (s - 3600h)/60);
             @sprintf("%02d:%02d:00", h, m))

"Escribe <name>.toml + <name>.pbs y devuelve (pbs=ruta, wall=segundos)."
function emit_job(name, cfg, outdir, ppn, base_dir, deltas, dt, est, wall)
    toml_path = joinpath(outdir, "$(name).toml")
    pbs_path  = joinpath(outdir, "$(name).pbs")
    open(io -> TOML.print(io, cfg), toml_path, "w")

    open(pbs_path, "w") do io
        println(io, """
        #!/bin/bash
        #PBS -N $(name)
        #PBS -l nodes=1:ppn=$(ppn)
        #PBS -l walltime=$(hhmmss(wall))
        #PBS -o results/logs/$(name)_\$PBS_JOBID.log
        #PBS -j oe

        cd \$PBS_O_WORKDIR
        mkdir -p results/logs "$(base_dir)"

        echo "================================"
        echo "$(name)  job \$PBS_JOBID  on \$(hostname)"
        echo "deltas: $(deltas)   ppn=$(ppn)   dt=$(dt)"
        echo "estimado ~$(round(est/3600; digits=1)) h, walltime $(hhmmss(wall))"
        echo "started \$(date '+%F %T')"
        echo "================================"

        julia --project=. -t $(ppn) scripts/compute_ensemble_Nsweep.jl "$(toml_path)"
        EXIT=\$?

        # El .jld2 real es la única prueba de éxito: un exit 0 sin fichero es un fallo.
        RESULT=\$(find "$(base_dir)" -name '*.jld2' 2>/dev/null | head -1)
        if [ \$EXIT -eq 0 ] && [ -z "\$RESULT" ]; then
          echo "✗ exit 0 pero NO se escribió .jld2"
          EXIT=1
        fi

        echo "finished \$(date '+%F %T')"
        if [ \$EXIT -eq 0 ]; then
          echo "✓ OK  -> \$RESULT"; ls -lh "\$RESULT"
          touch "results/logs/$(name)_\$PBS_JOBID.SUCCESS"
        else
          echo "✗ FALLO (exit \$EXIT)"
          touch "results/logs/$(name)_\$PBS_JOBID.FAILED"
        fi
        exit \$EXIT
        """)
    end
    chmod(pbs_path, 0o755)
    (pbs = pbs_path, wall = wall)
end

function main()
    isempty(ARGS) && error("Uso: julia scripts/generate_pbs_nsweep.jl <config.toml> [--ns-per-step-site=C] [--safety=F]")
    cfgpath = ARGS[1]
    C       = parse(Float64, getopt(ARGS, "ns-per-step-site", "67"))
    safety  = parse(Float64, getopt(ARGS, "safety", "2.5"))
    outdir    = getopt(ARGS, "outdir", "jobs")
    per_delta = "--per-delta" in ARGS
    maxwall   = 3600 * parse(Float64, getopt(ARGS, "max-walltime-h", "48"))
    mkpath(outdir)
    jobs = NamedTuple{(:pbs, :wall), Tuple{String, Float64}}[]

    d    = TOML.parsefile(cfgpath)
    phys = d["physics"]; sim = d["simulation"]
    haskey(phys, "N_values") || error("El config debe tener 'N_values' en [physics]")

    N_values = Vector{Int}(phys["N_values"])
    deltas   = Float64.(phys["delta_values"])
    boundary = Symbol(phys["boundary"])
    stype    = phys["system_type"]
    dt       = Float64(sim["DT"])
    base_dir = d["output"]["base_dir"]
    stmax    = haskey(sim, "scaled_t_max") ? Float64(sim["scaled_t_max"]) : nothing

    @printf("=== Generador PBS nsweep ===\nconfig: %s\nC = %.1f ns/paso/sitio, safety = %.1f, dt = %g\n\n",
            cfgpath, C, safety, dt)
    isnothing(stmax) && println("Nota: sin 'scaled_t_max'; se usa TMAX fijo del config.\n")

    @printf("%6s %5s %12s %12s %12s %12s\n", "N", "ppn", "TMAX(max)", "t_caso", "estimado", "walltime")
    total = 0.0
    for N in N_values
        # El Δκ más caro manda: ω₂ menor ⇒ TMAX mayor.
        tmaxs = [isnothing(stmax) ? Float64(sim["TMAX"]) :
                 stmax * 2π / omega_ref(N, dk, boundary, stype) for dk in deltas]
        TMAX_max = maximum(tmaxs)

        t_caso = C * 1e-9 * N * (TMAX_max / dt)

        # Un job por (N,Δκ) con 1 core evita por completo la contención entre hilos:
        # el walltime baja de t_caso·1.7·safety a t_caso·safety. Es lo que hace que
        # N=256 quepa por debajo del límite típico de 24 h de las colas.
        if per_delta
            for (j, dk) in enumerate(deltas)
                tmax_j = tmaxs[j]
                tj     = C * 1e-9 * N * (tmax_j / dt)
                wj     = tj * safety
                total += tj

                cfgj = deepcopy(d)
                delete!(cfgj["physics"], "N_values")
                cfgj["physics"]["N_values"]     = [N]
                cfgj["physics"]["N"]            = N
                cfgj["physics"]["delta_values"] = [dk]
                tag = replace(string(dk), "." => "p")
                cfgj["output"]["base_dir"] = "$(base_dir)_N$(N)_d$(tag)"
                push!(jobs, emit_job("nsweep_N$(N)_d$(tag)", cfgj, outdir, 1,
                                     cfgj["output"]["base_dir"], [dk], dt, tj, wj))
                @printf("%6d %5d %12.3e %12s %12s %12s   Δκ=%.2f\n",
                        N, 1, tmax_j, hhmmss(tj), hhmmss(tj), hhmmss(wj), dk)
            end
            continue
        end

        ppn    = min(length(deltas), 8)
        rondas = ceil(Int, length(deltas) / ppn)
        est    = t_caso * rondas * CONTENTION
        wall   = est * safety
        total += est

        cfg = deepcopy(d)
        delete!(cfg["physics"], "N_values")
        cfg["physics"]["N"] = N              # informativo; el script usa N_values
        cfg["physics"]["N_values"] = [N]
        cfg["output"]["base_dir"] = "$(base_dir)_N$(N)"

        push!(jobs, emit_job("nsweep_N$(N)", cfg, outdir, ppn,
                             cfg["output"]["base_dir"], deltas, dt, est, wall))

        @printf("%6d %5d %12.3e %12s %12s %12s\n", N, ppn, TMAX_max,
                hhmmss(t_caso), hhmmss(est), hhmmss(wall))
    end

    @printf("\nTotal si se lanzan en serie: ~%.1f h.  En cola paralela: el job más largo.\n", total/3600)
    if maxwall > 0 && any(j -> j.wall > maxwall, jobs)
        @printf("\n⚠  %d job(s) superan el límite de %g h de la cola.\n",
                count(j -> j.wall > maxwall, jobs), maxwall/3600)
        println("   Relanzar con --per-delta: 1 core por (N,Δκ) elimina la contención")
        println("   (×1.7) y baja el walltime lo suficiente para entrar.")
    end
    println("\nLanzar:")
    for j in jobs
        println("  qsub $(j.pbs)")
    end
    println("\n⚠  C=$(C) ns es del portátil de desarrollo. Medir en el clúster primero:")
    println("   julia --project=. scripts/calibrate_cost.jl")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
