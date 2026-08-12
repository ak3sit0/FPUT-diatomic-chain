"""
    calibrate_cost.jl

Mide la constante de coste de esta máquina para el integrador FPUT:

    wall_por_trayectoria ≈ C · N · (TMAX / dt)

C se expresa en ns/paso/sitio y es la ÚNICA cantidad dependiente de la máquina.
Medida en el portátil de desarrollo: C ≈ 67 ns (tras la reescritura del kernel).
En el clúster será distinta — ejecutar esto ahí ANTES de dimensionar los walltime.

Uso:
  julia --project=. scripts/calibrate_cost.jl            # N por defecto
  julia --project=. scripts/calibrate_cost.jl 64,128,256

Notas de método (importan, se midieron mal la primera vez):
  - Hay que CALENTAR por cada (N, dt): una sonda en frío mide la compilación JIT,
    no el cálculo. Una sonda a T=500 llegó a marcar 301 µs/paso frente a 30 reales.
  - El coste por paso es independiente de T (verificado en un rango de 64×), así que
    una sonda corta (T≈2000) predice exactamente una tirada de producción larga.
  - BLAS a 1 hilo: si no, cada hilo de Julia lanza los suyos y contaminan la medida.
"""

using LinearAlgebra, Printf, Statistics, Random
include("../src/fput_core.jl");        using .FPUTCore
include("../src/fput_fast_runner.jl"); using .FPUTFastRunner

const T_PROBE = 2000.0      # ≥2000: por debajo, la compilación domina
const N_REPS  = 3

function probe(N::Int, dt::Float64)
    sp  = SystemParams(N, 0.1, 0.0, 0.1, 0.0, :periodic)
    rng = Xoshiro(1)
    q0  = 0.01 .* randn(rng, N)
    v0  = zeros(N)
    saveat = 0.0:(dt * 3N):T_PROBE

    solve_fput(sp, q0, v0, (0.0, 50.0), dt; saveat=[0.0, 50.0])       # compilar
    solve_fput(sp, q0, v0, (0.0, T_PROBE), dt; saveat=saveat)         # calentar cachés

    ts = Float64[]
    for _ in 1:N_REPS
        GC.gc()
        push!(ts, @elapsed solve_fput(sp, q0, v0, (0.0, T_PROBE), dt; saveat=saveat))
    end
    nsteps = round(Int, T_PROBE / dt)
    (min = 1e9 * minimum(ts) / nsteps / N,
     med = 1e9 * median(ts)  / nsteps / N,
     spread = 100 * (maximum(ts) - minimum(ts)) / minimum(ts))
end

function main()
    BLAS.set_num_threads(1)
    N_values = isempty(ARGS) ? [64, 128, 256, 512] :
               parse.(Int, split(ARGS[1], ","))

    println("=== Calibración de coste (host: $(gethostname())) ===")
    println("Julia $(VERSION), BLAS threads = 1, sonda T=$(T_PROBE), $(N_REPS) repeticiones\n")
    @printf("%6s %14s %12s %9s\n", "N", "ns/paso/sitio", "mediana", "dispersión")

    cs = Float64[]
    for N in N_values
        r = probe(N, 0.2)
        push!(cs, r.min)
        @printf("%6d %14.1f %12.1f %8.0f%%\n", N, r.min, r.med, r.spread)
    end

    # Para dimensionar walltime se usa el MÁXIMO, no la mediana: sobrestimar la duración
    # es inocuo (se pide de más), subestimarla mata el job al llegar al límite de la cola.
    C   = maximum(cs)
    rel = (maximum(cs) - minimum(cs)) / minimum(cs)
    println()
    @printf("C (máximo, para walltime) = %.1f ns/paso/sitio\n", C)
    @printf("mediana = %.1f, rango = %.0f%%\n", median(cs), 100rel)

    if rel > 0.25
        # El SIGNO de la tendencia importa: creciente = pared de memoria (grave);
        # decreciente = overhead fijo por paso amortizándose (benigno).
        if cs[end] > cs[1]
            println("\n⚠  ns/paso/sitio CRECE con N: el modelo lineal se rompe (probable")
            println("   límite de memoria/caché). NO extrapolar a N mayores; medir cada N.")
        else
            println("\n✓  ns/paso/sitio DECRECE con N: es el overhead fijo por paso (las 18")
            println("   etapas de KahanLi8) amortizándose. Benigno — el coste por sitio sólo")
            println("   mejora al crecer N, así que usar el máximo va sobrado.")
        end
    else
        println("\n✓  Plano en N ⇒ coste = C·N·TMAX/dt es fiable para extrapolar.")
    end

    mkpath("results/logs")
    write("results/logs/cost_constant.txt", string(round(C; digits=1)))
    println("\nGuardado en results/logs/cost_constant.txt")
    println("Pasar a generate_pbs_nsweep.jl con:  --ns-per-step-site=$(round(C; digits=1))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
