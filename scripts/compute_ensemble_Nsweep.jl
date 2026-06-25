"""
    compute_ensemble_Nsweep.jl

Barrido sobre múltiples N con parámetros fijos (α, Δκ).
Un hilo por N, una realización por N (init_type = "mode" o "band_ensemble").
Salida: un solo JLD2 con todas las N (retrocompatible con plot_entropy_N_timeseries.jl).

Usage:
  julia --project=. -t 4 scripts/compute_ensemble_Nsweep.jl configs/cases/nsweep_test_quick.toml
"""

using TOML, JLD2, Dates, Base.Threads, Statistics, LinearAlgebra, Random
include("../src/fput_core.jl");        using .FPUTCore
include("../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../src/fput_analysis.jl");    using .FPUTAnalysis

# ── Funciones copiadas de compute_ensemble.jl ────────────────────────────

function derive_band_indices(branch::String, N::Int, freq::Vector{Float64},
                              boundary::Symbol;
                              k_band_start::Union{Int,Nothing}=nothing,
                              k_band_end::Union{Int,Nothing}=nothing)
    if !isnothing(k_band_start) && !isnothing(k_band_end)
        return collect(k_band_start:k_band_end)
    end

    if boundary == :periodic
        diffs = diff(freq)
        gap_idx = argmax(diffs[2:end]) + 1
        if branch == "acoustic"
            return collect(2:gap_idx)
        elseif branch == "optical"
            return collect(gap_idx+1:N)
        else
            error("branch debe ser 'acoustic' u 'optical', recibido: '$branch'")
        end
    else  # :fixed
        branch == "acoustic" && return collect(1:N)
        error("Frontera fija no tiene rama óptica distinguible")
    end
end

function band_phase_ic(k_band::AbstractVector{Int}, E_total::Float64, N::Int,
                       freq::Vector{Float64}, V::Matrix{Float64},
                       m::Vector{Float64}, seed::Int)
    rng    = Random.Xoshiro(seed)
    E_per  = E_total / length(k_band)
    Q      = zeros(N)
    P      = zeros(N)

    for j in k_band
        freq[j] < 1e-10 && continue
        φ    = rand(rng) * 2π
        A    = sqrt(2 * E_per)
        Q[j] =  A / freq[j] * cos(φ)
        P[j] = -A            * sin(φ)
    end

    inv_sqrt_m = 1.0 ./ sqrt.(m)
    q0 = (V * Q) .* inv_sqrt_m
    v0 = (V * P) .* inv_sqrt_m

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

function run_single_realization(sp, q0, v0, freq, V, m, target_mode_idx, k_ac, k_opt, cfg, label, N::Int)
    q_cur = copy(q0)
    v_cur = copy(v0)
    T_total        = Float64[]
    modal_E_blocks = Vector{Matrix{Float64}}()
    E_ac_blocks    = Vector{Vector{Float64}}()
    E_opt_blocks   = Vector{Vector{Float64}}()
    t_cur          = 0.0

    while t_cur < cfg.TMAX
        t_next = min(t_cur + cfg.T_block, cfg.TMAX)
        saveat = t_cur:cfg.save_every*cfg.DT:t_next

        Qb, Vb, Tb, _, _ = FPUTFastRunner.solve_fput(sp, q_cur, v_cur,
                                                       (t_cur, t_next), cfg.DT;
                                                       saveat=saveat)

        Tb = Tb .* freq[target_mode_idx] ./ (2π)

        Qmat = size(Qb,1) == N ? Float64.(Qb) : Float64.(Qb')
        Vmat = size(Vb,1) == N ? Float64.(Vb) : Float64.(Vb')

        modal_Eb = FPUTAnalysis.compute_modal_energies(Qmat, Vmat, freq, V, m)

        E_cur = sum(modal_Eb[:, end])
        if E_cur > 100 * cfg.E_total || isnan(E_cur) || isinf(E_cur)
            println("  $label INESTABILIDAD detectada (E=$(round(E_cur; sigdigits=3)) >> E_total=$(cfg.E_total)). Abortando.")
            return nothing
        end

        E_ac_b  = vec(sum(modal_Eb[k_ac,  :], dims=1))
        E_opt_b = vec(sum(modal_Eb[k_opt, :], dims=1))

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

# ── Configuración ──────────────────────────────────────────────────────

function build_nsweep_config(path::String)
    d    = TOML.parsefile(path)
    phys = d["physics"]
    sim  = d["simulation"]
    out  = d["output"]

    # Leer N_values
    N_values = if haskey(phys, "N_values")
        Vector{Int}(phys["N_values"])
    else
        error("Config debe tener 'N_values' en [physics]")
    end

    E_total = if haskey(phys, "energy_density")
        Float64(phys["energy_density"]) * first(N_values)  # Usar primer N para calcular energía
    else
        Float64(get(phys, "initial_energy", 0.45))
    end

    return (
        N_values      = N_values,
        boundary      = Symbol(phys["boundary"]),
        system_type   = phys["system_type"],
        nonlinear     = Symbol(phys["nonlinear"]),
        param_values  = Float64.(phys["param_values"]),
        delta_values  = Float64.(phys["delta_values"]),
        E_total       = E_total,
        init_type     = get(phys, "init_type", "mode"),
        branch        = get(phys, "branch", "acoustic"),
        n_real        = Int(get(phys, "n_real", 1)),
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

# ── Función principal: correr un N ─────────────────────────────────────

function run_for_N(N::Int, pval::Float64, delta::Float64, cfg)
    """Ejecuta una realización para un N específico."""
    label = "N=$N p=$pval Δ=$delta"
    println("[$(label)]")

    # Ajustar energía para este N específico
    E_total = if haskey(cfg, :energy_density)
        get(cfg, :energy_density, 0.445) * N
    else
        cfg.E_total
    end

    sp = FPUTCore.SystemParams(
        N,
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

    # Seleccionar índices de banda
    if cfg.init_type == "band_ensemble"
        k_ac = derive_band_indices(cfg.branch, N, freq, cfg.boundary;
                                   k_band_start=cfg.k_band_start,
                                   k_band_end=cfg.k_band_end)
        other_branch = cfg.branch == "acoustic" ? "optical" : "acoustic"
        k_opt        = derive_band_indices(other_branch, N, freq, cfg.boundary)
        k_band       = k_ac
        target_mode_idx = k_ac[1]
    else  # "mode"
        target_mode_idx = cfg.boundary == :fixed ? 1 : 2
        k_band = [target_mode_idx]
        k_ac   = collect(1:div(N, 2))
        k_opt  = collect(div(N, 2)+1:N)
    end

    # Condición inicial
    if cfg.init_type == "band_ensemble"
        seed = cfg.seed_base + 1
        q0, v0 = band_phase_ic(k_band, E_total, N, freq, V, m, seed)
    else  # "mode"
        U         = Diagonal(1.0 ./ sqrt.(m)) * V
        amplitude = sqrt(2 * E_total) / freq[target_mode_idx]
        q0        = amplitude .* U[:, target_mode_idx]
        v0        = zeros(N)
    end

    # Ejecutar
    result = run_single_realization(sp, q0, v0, freq, V, m,
                                    target_mode_idx, k_ac, k_opt, cfg, label, N)

    if isnothing(result)
        println("  ✗ Fallo")
        return nothing
    end

    println("  ✓ Completado: S_final=$(round(result.entropy[end]; digits=4))")

    # Retornar entrada compatible con plot_entropy_N_timeseries.jl
    return (
        N         = N,
        param     = pval,
        Delta     = delta,
        scaled_t  = result.scaled_t,
        modal_E   = result.modal_E,
        entropy   = result.entropy,
        E_total   = E_total,
        init_type = cfg.init_type,
        branch    = cfg.branch,
    )
end

# ── Main ───────────────────────────────────────────────────────────────

function main()
    isempty(ARGS) && error("Usage: julia scripts/compute_ensemble_Nsweep.jl <config.toml>")
    cfg = build_nsweep_config(ARGS[1])
    mkpath(cfg.base_dir)

    # Parámetros fijos (debe haber solo 1 param y 1 delta)
    pval  = cfg.param_values[1]
    delta = cfg.delta_values[1]

    println("=== Nsweep ===")
    println("Parámetros: α=$pval, Δκ=$delta")
    println("N_values: $(cfg.N_values)")
    println("init_type: $(cfg.init_type)")
    println()

    # Paralelismo: 1 hilo por N
    N_vals = cfg.N_values
    results = Vector{Any}(undef, length(N_vals))

    @threads for i in eachindex(N_vals)
        N = N_vals[i]
        results[i] = try
            run_for_N(N, pval, delta, cfg)
        catch e
            println("[N=$N] ERROR: $e\n$(sprint(showerror, e, catch_backtrace()))")
            nothing
        end
    end

    valid = filter(!isnothing, results)

    outpath = joinpath(cfg.base_dir, "nsweep_results_$(Dates.today()).jld2")
    jldsave(outpath; results=valid, config=ARGS[1])
    println("\nGuardado: $outpath  ($(length(valid))/$(length(N_vals)) casos válidos)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
