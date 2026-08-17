"""
    compute_ensemble.jl

Barrido en delta_values con ensamble de fases por banda selectiva.
Hereda el patrón de bloques temporales de compute_trajectories.jl.

Usage:
  julia --project=. scripts/compute_ensemble.jl configs/production/ensemble.toml

Convención de energía:
  Si el TOML tiene `energy_density`, E_total = N * energy_density  (nueva).
  Si solo tiene `initial_energy`,   E_total = initial_energy        (legado).
"""

using TOML, JLD2, Dates, Base.Threads, Statistics, LinearAlgebra, Random
include("../../src/fput_core.jl");        using .FPUTCore
include("../../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../../src/fput_analysis.jl");    using .FPUTAnalysis
include("../../src/block_integration.jl"); using .BlockIntegration

# derive_band_indices y band_phase_ic viven en FPUTCore (src/fput_core.jl)

# ── Configuración ─────────────────────────────────────────────────────────────

function build_ensemble_config(path::String)
    d    = TOML.parsefile(path)
    phys = d["physics"]
    sim  = d["simulation"]
    out  = d["output"]

    N = Int(phys["N"])

    # Convención de energía: energy_density (nueva) tiene prioridad sobre initial_energy (legado)
    E_total = if haskey(phys, "energy_density")
        Float64(phys["energy_density"]) * N
    else
        Float64(get(phys, "initial_energy", 0.45))
    end

    return (
        N             = N,
        boundary      = Symbol(phys["boundary"]),
        system_type   = phys["system_type"],
        nonlinear     = Symbol(phys["nonlinear"]),
        param_values  = Float64.(phys["param_values"]),
        delta_values  = Float64.(phys["delta_values"]),
        E_total       = E_total,
        init_type     = get(phys, "init_type", "mode"),
        branch        = get(phys, "branch", "acoustic"),
        n_real        = Int(get(phys, "n_real", 8)),
        seed_base     = Int(get(phys, "seed_base", 42)),
        k_band_start  = haskey(phys, "k_band_start") ? Int(phys["k_band_start"]) : nothing,
        k_band_end    = haskey(phys, "k_band_end")   ? Int(phys["k_band_end"])   : nothing,
        TMAX          = Float64(sim["TMAX"]),
        T_block       = Float64(sim["T_block"]),
        DT            = Float64(sim["DT"]),
        save_every    = Int(sim["save_every"]),
        downsample    = Int(sim["downsample"]),
        entropy_delta = Float64(get(sim, "entropy_delta", 0.6)),
        debug         = Bool(get(sim, "debug", false)),
        base_dir      = out["base_dir"],
    )
end

# ── Una realización (integración en bloques) ──────────────────────────────────

# k_ac, k_opt: índices de modos acústicos y ópticos (para diagnóstico interbanda)
function run_single_realization(sp, q0, v0, freq, V, m, target_mode_idx, k_ac, k_opt, cfg, label)
    E_ac_blocks  = Vector{Vector{Float64}}()   # energía total banda acústica vs t
    E_opt_blocks = Vector{Vector{Float64}}()   # energía total banda óptica vs t

    on_block = (modal_Eb, idx_ds) -> begin
        push!(E_ac_blocks,  vec(sum(modal_Eb[k_ac,  idx_ds], dims=1)))
        push!(E_opt_blocks, vec(sum(modal_Eb[k_opt, idx_ds], dims=1)))
    end
    # Detección de inestabilidad: la energía total no debe alejarse >100× del valor inicial
    check_instability = modal_Eb -> begin
        E_cur = sum(modal_Eb[:, end])
        unstable = E_cur > 100 * cfg.E_total || isnan(E_cur) || isinf(E_cur)
        unstable && println("  $label INESTABILIDAD detectada (E=$(round(E_cur; sigdigits=3)) >> E_total=$(cfg.E_total)). Abortando.")
        unstable
    end

    result = integrate_in_blocks(sp, q0, v0, cfg.N, freq, V, m, target_mode_idx;
                                  TMAX=cfg.TMAX, T_block=cfg.T_block, DT=cfg.DT,
                                  save_every=cfg.save_every, downsample=cfg.downsample,
                                  debug=cfg.debug, label=label,
                                  on_block=on_block, check_instability=check_instability)

    result.aborted && return nothing

    entropy = FPUTAnalysis.spectral_entropy(result.modal_E, cfg.entropy_delta)
    E_ac    = isempty(E_ac_blocks)  ? Float64[] : reduce(vcat, E_ac_blocks)
    E_opt   = isempty(E_opt_blocks) ? Float64[] : reduce(vcat, E_opt_blocks)
    return (scaled_t=result.scaled_t, modal_E=result.modal_E, entropy=entropy,
            E_acoustic=E_ac, E_optical=E_opt)
end

# ── Tiempo de termalización (Opción A: umbral por realización) ────────────────

"""
    compute_T_therm(entropy_realizations, scaled_t, k_band, N) -> NamedTuple

Estima T_therm por realización via umbral en la entropía normalizada,
luego calcula estadísticas del ensamble. Umbral f = 1 - 1/e (tiempo de escala
natural: equivale al tiempo de relajación para crecimiento exponencial puro).

S0 = log(|k_band|)  entropía inicial teórica (banda uniforme)
Seq = log(N)         equipartición total
"""
function compute_T_therm(entropy_realizations, scaled_t, k_band, N)
    f  = 1 - 1/ℯ
    S0  = log(length(k_band))
    Seq = log(N)
    ΔS  = Seq - S0

    T_vec = map(entropy_realizations) do S_r
        nt     = min(length(S_r), length(scaled_t))
        S_norm = (S_r[1:nt] .- S0) ./ ΔS
        idx    = findfirst(S_norm .> f)
        isnothing(idx) ? Inf : scaled_t[idx]
    end

    finitos = filter(isfinite, T_vec)
    n_fin   = length(finitos)
    return (
        T_therm_mean   = isempty(finitos) ? Inf : mean(finitos),
        T_therm_std    = (n_fin > 1)      ? std(finitos) : NaN,
        T_therm_median = isempty(finitos) ? Inf : median(finitos),
        T_therm_vec    = T_vec,
        frac_therm     = n_fin / length(T_vec),
        n_therm        = n_fin,
        threshold_f    = f,
    )
end

# ── Un caso completo (ensamble) ───────────────────────────────────────────────

function run_ensemble_case(case_idx, pval, delta, cfg)
    println("[case $case_idx] p=$pval  Δ=$delta  (init=$(cfg.init_type), n_real=$(cfg.n_real))")

    sp = FPUTCore.SystemParams(
        cfg.N,
        cfg.system_type == "springs" ? delta : 0.0,
        cfg.system_type == "masses"  ? delta : 0.0,
        cfg.nonlinear == :alpha ? pval : 0.0,
        cfg.nonlinear == :beta  ? pval : 0.0,
        cfg.boundary
    )

    k, m    = FPUTCore.make_system(sp)
    freq, V = FPUTCore.find_normal_modes(k, m, cfg.boundary)
    idx_s   = sortperm(freq)
    freq    = freq[idx_s]
    V       = V[:, idx_s]

    # Selección de banda y modo de referencia para la escala temporal
    if cfg.init_type == "band_ensemble"
        k_ac = derive_band_indices(cfg.branch, cfg.N, freq, cfg.boundary;
                                   k_band_start=cfg.k_band_start,
                                   k_band_end=cfg.k_band_end)
        other_branch = cfg.branch == "acoustic" ? "optical" : "acoustic"
        k_opt        = derive_band_indices(other_branch, cfg.N, freq, cfg.boundary)
        k_band       = k_ac
        target_mode_idx = k_ac[1]
        println("  Banda '$(cfg.branch)': modos $(k_ac[1])..$(k_ac[end]) ($(length(k_ac)) modos)")
        println("  Banda complementaria: modos $(k_opt[1])..$(k_opt[end]) ($(length(k_opt)) modos)")
    else  # "mode" — backward compat
        target_mode_idx = cfg.boundary == :fixed ? 1 : 2
        k_band = [target_mode_idx]
        k_ac   = collect(1:div(cfg.N, 2))
        k_opt  = collect(div(cfg.N, 2)+1:cfg.N)
    end

    n_real = cfg.n_real

    # Acumuladores
    entropy_realizations  = Vector{Vector{Float64}}()
    E_optical_realizations = Vector{Vector{Float64}}()  # nueva: guardar E_opt de cada realización
    modal_E_mean_accum    = nothing
    E_ac_mean_accum       = nothing
    E_opt_mean_accum      = nothing
    T_ref                 = nothing
    n_ok                  = 0   # realizaciones estables

    # Las n_real realizaciones son independientes (sólo leen sp/freq/V/m; cada
    # solve_fput reserva su propio buffer), así que se lanzan en paralelo. Antes el
    # bucle era secuencial y el único paralelismo estaba en (param, Δκ) — 9 tareas —,
    # así que con ppn=16 sobraban 7 cores y pedir un nodo mayor no servía de nada.
    # Medido: 8 realizaciones a N=256, 11.4 s en serie → 2.2 s en 8 hilos (5.1×, GC 0%).
    # Se lanzan POR LOTES de nthreads() en vez de las n_real de golpe. Lanzarlas todas
    # deja vivos los n_real resultados a la vez hasta terminar la reducción: con
    # n_real=100 y modal_E de 13 MB son ~1.3 GB por caso. Por lotes el pico queda
    # acotado a nthreads() resultados (~260 MB con ppn=20), y no se pierde tiempo:
    # con n_real múltiplo de ppn los lotes coinciden con las rondas que ya haría.
    batch = max(nthreads(), 1)

    for batch_start in 1:batch:n_real
        rs    = batch_start:min(batch_start + batch - 1, n_real)
        tasks = Vector{Any}(undef, length(rs))   # Any: se ponen a nothing al liberar
        for (i, r) in enumerate(rs)
            tasks[i] = Threads.@spawn begin
                seed = cfg.seed_base + r

                q0, v0 = if cfg.init_type == "band_ensemble"
                    band_phase_ic(k_band, cfg.E_total, cfg.N, freq, V, m, seed)
                else
                    U         = Diagonal(1.0 ./ sqrt.(m)) * V
                    amplitude = sqrt(2 * cfg.E_total) / freq[target_mode_idx]
                    (amplitude .* U[:, target_mode_idx], zeros(cfg.N))
                end

                run_single_realization(sp, q0, v0, freq, V, m,
                                       target_mode_idx, k_ac, k_opt, cfg,
                                       "case=$case_idx r=$r/$n_real")
            end
        end

    # La reducción se hace en serie y en orden de r, para que el resultado NO dependa
    # del orden en que terminen los hilos (reproducibilidad bit a bit con las semillas).
    for (i, r) in enumerate(rs)
        result = fetch(tasks[i])
        tasks[i] = nothing        # liberar el resultado en cuanto se acumula

        # Descartar realizaciones inestables
        isnothing(result) && (println("  [r=$r] descartada (inestable)"); continue)

        n_ok += 1
        push!(entropy_realizations, result.entropy)
        push!(E_optical_realizations, result.E_optical)  # nueva línea
        T_ref = result.scaled_t

        # Acumulación incremental con peso 1/n_real (se renormaliza al final si n_ok < n_real)
        if isnothing(modal_E_mean_accum)
            modal_E_mean_accum = copy(result.modal_E)
            E_ac_mean_accum    = copy(result.E_acoustic)
            E_opt_mean_accum   = copy(result.E_optical)
        else
            nt = min(size(result.modal_E, 2), size(modal_E_mean_accum, 2))
            modal_E_mean_accum = modal_E_mean_accum[:, 1:nt] .+ result.modal_E[:, 1:nt]
            E_ac_mean_accum    = E_ac_mean_accum[1:nt]       .+ result.E_acoustic[1:nt]
            E_opt_mean_accum   = E_opt_mean_accum[1:nt]      .+ result.E_optical[1:nt]
        end

        println("  [r=$r] S_final=$(round(result.entropy[end]; digits=4))  " *
                "E_ac=$(round(result.E_acoustic[end]; digits=4))  " *
                "E_opt=$(round(result.E_optical[end]; digits=4))")
    end
    end  # fin del lote

    if n_ok == 0
        println("[case $case_idx] Todas las realizaciones inestables. Descartando caso.")
        return nothing
    end

    # Normalizar por número de realizaciones estables
    modal_E_mean_accum ./= n_ok
    E_ac_mean_accum    ./= n_ok
    E_opt_mean_accum   ./= n_ok
    n_ok < n_real && println("  Advertencia: solo $n_ok/$n_real realizaciones estables.")

    # Estadísticas del ensamble sobre realizaciones estables
    nt_min   = minimum(length(e) for e in entropy_realizations)
    ent_mat  = reduce(hcat, [e[1:nt_min] for e in entropy_realizations])'
    entropy_mean = vec(mean(ent_mat, dims=1))
    entropy_std  = vec(std(ent_mat,  dims=1))

    ttherm = compute_T_therm(entropy_realizations, T_ref, k_band, cfg.N)

    println("[case $case_idx] Finalizado. S̄=$(round(entropy_mean[end];digits=4)) ± $(round(entropy_std[end];digits=4))  " *
            "T_therm=$(round(ttherm.T_therm_mean; sigdigits=3)) ± $(round(ttherm.T_therm_std; sigdigits=2))  " *
            "frac_therm=$(ttherm.n_therm)/$(n_real)")

    return (
        param                  = pval,
        Delta                  = delta,
        scaled_t               = T_ref,
        entropy_mean           = entropy_mean,
        entropy_std            = entropy_std,
        modal_E_mean           = modal_E_mean_accum,
        E_acoustic_mean        = E_ac_mean_accum,
        E_optical_mean         = E_opt_mean_accum,
        entropy_realizations   = entropy_realizations,
        E_optical_realizations = E_optical_realizations,
        T_therm_mean           = ttherm.T_therm_mean,
        T_therm_std            = ttherm.T_therm_std,
        T_therm_median         = ttherm.T_therm_median,
        T_therm_vec            = ttherm.T_therm_vec,
        frac_therm             = ttherm.frac_therm,
        n_therm                = ttherm.n_therm,
        n_real                 = n_real,
        seed_base              = cfg.seed_base,
        branch                 = cfg.branch,
        k_band                 = collect(k_band),
        E_total                = cfg.E_total,
    )
end

# ── Main ──────────────────────────────────────────────────────────────────────

function main()
    isempty(ARGS) && error("Usage: julia scripts/compute_ensemble.jl <config.toml>")
    cfg = build_ensemble_config(ARGS[1])
    mkpath(cfg.base_dir)

    pairs   = collect(Iterators.product(cfg.param_values, cfg.delta_values))
    M       = length(pairs)
    results = Vector{Any}(undef, M)

    # Los casos van EN SERIE; el paralelismo está ahora dentro, sobre las n_real
    # realizaciones (ver run_ensemble_case). Anidar @threads aquí con los @spawn de
    # dentro infrautilizaría los hilos: el bucle externo los tendría bloqueados
    # esperando. Con n_real ≫ ppn el bucle interno ya satura el nodo por sí solo.
    println("$(M) casos en serie, $(cfg.n_real) realizaciones en paralelo sobre $(nthreads()) hilos\n")
    for i in 1:M
        pval, delta = pairs[i]
        results[i]  = try
            run_ensemble_case(i, pval, delta, cfg)
        catch e
            println("[case $i] ERROR: $e\n$(sprint(showerror, e, catch_backtrace()))")
            nothing
        end
    end

    valid   = filter(!isnothing, results)
    outpath = joinpath(cfg.base_dir, "ensemble_results_$(Dates.today()).jld2")
    jldsave(outpath; results=valid, config=ARGS[1])
    println("Guardado: $outpath  ($(length(valid))/$(M) casos válidos)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
