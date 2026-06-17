"""
    compute_ensemble.jl

Barrido en delta_values con ensamble de fases por banda selectiva.
Hereda el patrón de bloques temporales de compute_trajectories.jl.

Usage:
  julia --project=. scripts/compute_ensemble.jl configs/cases/ensemble_test_quick.toml

Convención de energía:
  Si el TOML tiene `energy_density`, E_total = N * energy_density  (nueva).
  Si solo tiene `initial_energy`,   E_total = initial_energy        (legado).
"""

using TOML, JLD2, Dates, Base.Threads, Statistics, LinearAlgebra, Random
include("../src/fput_core.jl");        using .FPUTCore
include("../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../src/fput_analysis.jl");    using .FPUTAnalysis

# ── Funciones de inicialización ───────────────────────────────────────────────

"""
    derive_band_indices(branch, N, freq, boundary; k_band_start, k_band_end)

Devuelve Vector{Int} (1-based) de los índices de modo para la rama solicitada.
`freq` debe estar ordenado en forma ascendente antes de llamar esta función.
Si k_band_start/k_band_end están presentes, los usa directamente (override manual).
"""
function derive_band_indices(branch::String, N::Int, freq::Vector{Float64},
                              boundary::Symbol;
                              k_band_start::Union{Int,Nothing}=nothing,
                              k_band_end::Union{Int,Nothing}=nothing)
    if !isnothing(k_band_start) && !isnothing(k_band_end)
        return collect(k_band_start:k_band_end)
    end

    if boundary == :periodic
        # Buscar el mayor salto de frecuencia (band gap), ignorando modo de Goldstone (idx 1)
        diffs = diff(freq)
        gap_idx = argmax(diffs[2:end]) + 1   # índice donde ocurre el gap (en freq)
        if branch == "acoustic"
            return collect(2:gap_idx)         # excluir modo de traslación (ω≈0)
        elseif branch == "optical"
            return collect(gap_idx+1:N)
        else
            error("branch debe ser 'acoustic' o 'optical', recibido: '$branch'")
        end
    else  # :fixed
        branch == "acoustic" && return collect(1:N)
        error("Frontera fija no tiene rama óptica distinguible")
    end
end

"""
    band_phase_ic(k_band, E_total, N, freq, V, m, seed) -> (q0, v0)

Condición inicial de banda selectiva con fases aleatorias uniformes.

Energía E_total distribuida uniformemente sobre los modos en k_band:
  E_j = E_total / length(k_band)  para j ∈ k_band
  Q_j =  sqrt(2·E_j) / ω_j · cos(φ_j),  P_j = -sqrt(2·E_j) · sin(φ_j)
  q = (1/√m) · V · Q,   v = (1/√m) · V · P

Validación post-construcción con @assert (precisión de máquina).
"""
function band_phase_ic(k_band::AbstractVector{Int}, E_total::Float64, N::Int,
                       freq::Vector{Float64}, V::Matrix{Float64},
                       m::Vector{Float64}, seed::Int)
    rng    = Random.Xoshiro(seed)
    E_per  = E_total / length(k_band)
    Q      = zeros(N)
    P      = zeros(N)

    for j in k_band
        freq[j] < 1e-10 && continue   # modo de Goldstone (traslación, ω≈0)
        φ    = rand(rng) * 2π
        A    = sqrt(2 * E_per)
        Q[j] =  A / freq[j] * cos(φ)
        P[j] = -A            * sin(φ)
    end

    inv_sqrt_m = 1.0 ./ sqrt.(m)
    q0 = (V * Q) .* inv_sqrt_m
    v0 = (V * P) .* inv_sqrt_m

    # Proyectar de vuelta para verificar
    x        = sqrt.(m) .* q0
    vx       = sqrt.(m) .* v0
    Q_check  = V' * x
    P_check  = V' * vx
    E_check  = 0.5 .* (P_check.^2 .+ (freq.^2) .* Q_check.^2)
    E_out    = sum(E_check[setdiff(1:N, k_band)])
    E_in     = sum(E_check[k_band])
    tol      = 1e-8 * E_total

    @assert E_out < tol        "Fuga de energía fuera de banda: E_out=$(E_out) (tol=$(tol))"
    @assert abs(E_in - E_total) < tol "Normalización incorrecta: E_in=$(E_in) vs E_total=$(E_total)"

    return q0, v0
end

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
    q_cur = copy(q0)
    v_cur = copy(v0)
    T_total        = Float64[]
    modal_E_blocks = Vector{Matrix{Float64}}()
    E_ac_blocks    = Vector{Vector{Float64}}()   # energía total banda acústica vs t
    E_opt_blocks   = Vector{Vector{Float64}}()   # energía total banda óptica vs t
    t_cur          = 0.0

    while t_cur < cfg.TMAX
        t_next = min(t_cur + cfg.T_block, cfg.TMAX)
        saveat = t_cur:cfg.save_every*cfg.DT:t_next

        Qb, Vb, Tb, _, _ = FPUTFastRunner.solve_fput(sp, q_cur, v_cur,
                                                       (t_cur, t_next), cfg.DT;
                                                       saveat=saveat)

        # Escalar tiempo a ciclos del modo de referencia (misma convención que compute_trajectories)
        Tb = Tb .* freq[target_mode_idx] ./ (2π)

        # Normalizar forma: garantizar N × nt
        Qmat = size(Qb,1) == cfg.N ? Float64.(Qb) : Float64.(Qb')
        Vmat = size(Vb,1) == cfg.N ? Float64.(Vb) : Float64.(Vb')

        modal_Eb = FPUTAnalysis.compute_modal_energies(Qmat, Vmat, freq, V, m)

        # Detección de inestabilidad: la energía total no debe alejarse >100× del valor inicial
        E_cur = sum(modal_Eb[:, end])
        if E_cur > 100 * cfg.E_total || isnan(E_cur) || isinf(E_cur)
            println("  $label INESTABILIDAD detectada (E=$(round(E_cur; sigdigits=3)) >> E_total=$(cfg.E_total)). Abortando.")
            return nothing
        end

        # Diagnóstico interbanda: suma de energía por banda en cada instante
        E_ac_b  = vec(sum(modal_Eb[k_ac,  :], dims=1))
        E_opt_b = vec(sum(modal_Eb[k_opt, :], dims=1))

        # Append evitando duplicado en frontera de bloque
        if !isempty(T_total) && !isempty(Tb) && isapprox(T_total[end], Tb[1]; atol=1e-12, rtol=0)
            if length(Tb) > 1
                idx_ds = 2:cfg.downsample:length(Tb)
                append!(T_total, Float64.(Tb[idx_ds]))
                blk = Float64.(modal_Eb[:, idx_ds])
                size(blk,2) > 0 && push!(modal_E_blocks, blk)
                push!(E_ac_blocks,  Float64.(E_ac_b[idx_ds]))
                push!(E_opt_blocks, Float64.(E_opt_b[idx_ds]))
            end
        else
            idx_ds = 1:cfg.downsample:length(Tb)
            append!(T_total, Float64.(Tb[idx_ds]))
            blk = Float64.(modal_Eb[:, idx_ds])
            size(blk,2) > 0 && push!(modal_E_blocks, blk)
            push!(E_ac_blocks,  Float64.(E_ac_b[idx_ds]))
            push!(E_opt_blocks, Float64.(E_opt_b[idx_ds]))
        end

        q_cur .= Qmat[:, end]
        v_cur .= Vmat[:, end]
        t_cur  = t_next
        cfg.debug && println("  $label t=$(round(t_cur; digits=2)) / $(cfg.TMAX)")
    end

    modal_E = isempty(modal_E_blocks) ? zeros(cfg.N, 0) : reduce(hcat, modal_E_blocks)
    entropy = FPUTAnalysis.spectral_entropy(modal_E, cfg.entropy_delta)
    E_ac    = isempty(E_ac_blocks)  ? Float64[] : reduce(vcat, E_ac_blocks)
    E_opt   = isempty(E_opt_blocks) ? Float64[] : reduce(vcat, E_opt_blocks)
    return (scaled_t=T_total, modal_E=modal_E, entropy=entropy, E_acoustic=E_ac, E_optical=E_opt)
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
    modal_E_mean_accum    = nothing
    E_ac_mean_accum       = nothing
    E_opt_mean_accum      = nothing
    T_ref                 = nothing
    n_ok                  = 0   # realizaciones estables

    for r in 1:n_real
        seed = cfg.seed_base + r

        if cfg.init_type == "band_ensemble"
            q0, v0 = band_phase_ic(k_band, cfg.E_total, cfg.N, freq, V, m, seed)
        else
            U         = Diagonal(1.0 ./ sqrt.(m)) * V
            amplitude = sqrt(2 * cfg.E_total) / freq[target_mode_idx]
            q0        = amplitude .* U[:, target_mode_idx]
            v0        = zeros(cfg.N)
        end

        label  = "case=$case_idx r=$r/$n_real"
        result = run_single_realization(sp, q0, v0, freq, V, m,
                                        target_mode_idx, k_ac, k_opt, cfg, label)

        # Descartar realizaciones inestables
        isnothing(result) && (println("  [r=$r] descartada (inestable)"); continue)

        n_ok += 1
        push!(entropy_realizations, result.entropy)
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

    println("[case $case_idx] Finalizado. S̄=$(round(entropy_mean[end];digits=4)) ± $(round(entropy_std[end];digits=4))")

    return (
        param                = pval,
        Delta                = delta,
        scaled_t             = T_ref,
        entropy_mean         = entropy_mean,
        entropy_std          = entropy_std,
        modal_E_mean         = modal_E_mean_accum,
        E_acoustic_mean      = E_ac_mean_accum,
        E_optical_mean       = E_opt_mean_accum,
        entropy_realizations = entropy_realizations,
        n_real               = n_real,
        seed_base            = cfg.seed_base,
        branch               = cfg.branch,
        k_band               = collect(k_band),
        E_total              = cfg.E_total,
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

    # Paralelismo sobre (param, delta); el ensamble interno es secuencial
    @threads for i in 1:M
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
