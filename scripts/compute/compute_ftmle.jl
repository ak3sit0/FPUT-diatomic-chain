"""
    compute_ftmle.jl

Calcula el ftMLE vía el método de dos trayectorias de Benettin et al.
Soporta múltiples Δκ en paralelo (Threads.@threads).

Usage:
  julia -t <N> --project=. scripts/compute_ftmle.jl [config.toml] [T_max] [T_renorm] [delta1,delta2,...]

  delta1,delta2,...  lista separada por comas; si se omite usa todos los delta_values del TOML
  T_max              tiempo físico máximo (default 1e7)
  T_renorm           intervalo de renormalización (default 200.0)
"""

using TOML, JLD2, LinearAlgebra, Random, Dates

include("../../src/fput_core.jl");        using .FPUTCore
include("../../src/fput_fast_runner.jl"); using .FPUTFastRunner

# ── Lectura del config ──────────────────────────────────────────────────────

function read_config(path::String)
    d    = TOML.parsefile(path)
    phys = d["physics"]
    sim  = d["simulation"]
    out  = d["output"]
    return (
        N             = Int(phys["N"]),
        boundary      = Symbol(phys["boundary"]),
        system_type   = String(phys["system_type"]),
        nonlinear     = Symbol(phys["nonlinear"]),
        param_values  = Float64.(phys["param_values"]),
        delta_values  = Float64.(phys["delta_values"]),
        initial_energy = Float64(get(phys, "initial_energy", 0.445)),
        init_mode     = Int(get(phys, "init_mode", 2)),
        DT            = Float64(get(sim, "DT", 0.05)),
        base_dir      = out["base_dir"],
    )
end

# ── Construcción de la CI (modo único, determinista) ──────────────────────

"""
    mode_ic(N, delta_k, alpha, boundary, mode_idx, E) -> (q0, v0, freq, V, m, sp)

Condición inicial de modo puro: toda la energía en el modo `mode_idx`.
Replica exactamente la lógica de compute_trajectories.jl.
"""
function mode_ic(N, delta_k, alpha, boundary, mode_idx, E)
    sp   = SystemParams(N, delta_k, 0.0, alpha, 0.0, boundary)
    k, m = make_system(sp)
    freq, V = find_normal_modes(k, m, boundary)

    idx_sort = sortperm(freq)
    freq = freq[idx_sort]
    V    = V[:, idx_sort]

    amplitude = sqrt(2 * E) / freq[mode_idx]
    U  = Diagonal(1.0 ./ sqrt.(m)) * V
    q0 = amplitude .* U[:, mode_idx]
    v0 = zeros(N)
    return q0, v0, freq, V, m, sp
end

# ── Protocolo de Benettin ────────────────────────────────────────────────────

"""
    run_ftmle(cfg, delta_k; T_max, T_renorm, delta0, rng_seed) -> (t_vec, lambda_vec, omega_ref)

Integra dos trayectorias en chunks de T_renorm, renormalizando la separación
a delta0 en cada paso. Devuelve el vector de tiempos (físico) y λ(t) acumulado.
"""
function run_ftmle(cfg, delta_k;
                   T_max    = 80000.0,
                   T_renorm = 50.0,
                   delta0   = 1e-8,
                   rng_seed = 1234)

    alpha     = cfg.param_values[1]
    N         = cfg.N
    dt        = cfg.DT
    boundary  = cfg.boundary
    mode_idx  = cfg.init_mode
    E         = cfg.initial_energy

    q_ref, v_ref, freq, V, m, sp = mode_ic(N, delta_k, alpha, boundary, mode_idx, E)

    # Perturbación aleatoria normalizada en espacio de fase 2N
    rng = MersenneTwister(rng_seed) # Semilla fija para reproducibilidad de números aleatorios
    dv  = randn(rng, 2N) # Perturbación aleatoria normalizada en espacio de fase 2N
    dv ./= norm(dv) # Normalizar a 1
    q_per = q_ref .+ delta0 .* dv[1:N] # Perturbación inicial en q
    v_per = v_ref .+ delta0 .* dv[N+1:end] # Perturbación inicial en v

    n_chunks   = Int(ceil(T_max / T_renorm))
    t_vec      = Vector{Float64}(undef, n_chunks)
    lambda_vec = Vector{Float64}(undef, n_chunks)
    Lambda     = 0.0
    t          = 0.0

    println("ftMLE: N=$N  Δκ=$delta_k  α=$alpha  modo=$mode_idx  T_max=$T_max  chunks=$n_chunks")
    flush(stdout)

    for i in 1:n_chunks
        t_next = t + T_renorm

        # Integrar ambas trayectorias (guardar solo el estado final)
        Q1, V1, _, _, _ = solve_fput(sp, q_ref, v_ref, (t, t_next), dt; saveat=[t_next])
        Q2, V2, _, _, _ = solve_fput(sp, q_per, v_per, (t, t_next), dt; saveat=[t_next])

        q_ref = Q1[:, end]
        v_ref = V1[:, end]
        q2f   = Q2[:, end]
        v2f   = V2[:, end]

        # Separación en espacio de fase
        Δq = q2f .- q_ref
        Δv = v2f .- v_ref
        d  = norm(vcat(Δq, Δv))

        Lambda += log(d / delta0)
        t       = t_next
        lambda_t = Lambda / t

        t_vec[i]      = t
        lambda_vec[i] = lambda_t

        # Renormalizar: acercar la trayectoria perturbada de vuelta a delta0
        q_per = q_ref .+ (delta0 / d) .* Δq
        v_per = v_ref .+ (delta0 / d) .* Δv

        if mod(i, max(1, n_chunks ÷ 10)) == 0
            println("  chunk $i/$n_chunks  t=$(round(t; digits=1))  λ(t)=$(round(lambda_t; sigdigits=4))")
            flush(stdout)
        end
    end

    omega_ref = freq[mode_idx]
    return t_vec, lambda_vec, omega_ref
end

# ── Main ────────────────────────────────────────────────────────────────────

function save_result(cfg, toml_path, delta_k, T_max, t_vec, lambda_vec, omega_ref)
    t_cycles = t_vec .* omega_ref ./ (2π)
    outdir   = joinpath(cfg.base_dir, "ftmle")
    mkpath(outdir)
    bc_tag  = string(cfg.boundary)
    tag     = replace(string(delta_k), "." => "p")
    outpath = joinpath(outdir, "ftmle_$(bc_tag)_delta$(tag)_$(Dates.today()).jld2")
    jldsave(outpath;
        t_physical  = t_vec,
        t_cycles    = t_cycles,
        lambda      = lambda_vec,
        omega_ref   = omega_ref,
        delta_k     = delta_k,
        boundary    = string(cfg.boundary),
        init_mode   = cfg.init_mode,
        T_max       = T_max,
        config_path = toml_path,
    )
    println("Guardado: $outpath")
end

function main()
    toml_path = get(ARGS, 1, "configs/cases/periodic_N64_production.toml")
    T_max     = parse(Float64, get(ARGS, 2, "10000000.0"))
    T_renorm  = parse(Float64, get(ARGS, 3, "200.0"))

    isfile(toml_path) || error("Config not found: $toml_path")
    cfg = read_config(toml_path)

    # Arg 4: lista de deltas separada por comas, o todos los del TOML
    deltas = if length(ARGS) >= 4
        requested = parse.(Float64, split(ARGS[4], ","))
        for d in requested
            any(x -> isapprox(x, d; atol=1e-12), cfg.delta_values) ||
                error("delta=$d no está en cfg.delta_values=$(cfg.delta_values)")
        end
        requested
    else
        cfg.delta_values
    end

    println("Casos a calcular: Δκ = $deltas  ($(Threads.nthreads()) thread(s))")
    flush(stdout)

    Threads.@threads for delta_k in deltas
        t_vec, lambda_vec, omega_ref = run_ftmle(cfg, delta_k; T_max=T_max, T_renorm=T_renorm)
        save_result(cfg, toml_path, delta_k, T_max, t_vec, lambda_vec, omega_ref)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
