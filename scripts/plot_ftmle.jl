"""
    plot_ftmle.jl

Figura de tres paneles con eje x compartido (log10, en ciclos):
  - Superior:  S̄(t) del backup de entropía
  - Medio:     λ(t) (ftMLE) + referencia c·ln(t)/t  (log-log)
  - Inferior:  R(t) = λ(t) / [c·ln(t)/t]             (semi-log)

Usage:
  julia --project=. scripts/plot_ftmle.jl <ftmle.jld2> [entropy.jld2]
"""

using JLD2, Plots, LaTeXStrings, Statistics

function apply_global_plot_style!()
    default(titlefont = font(16), guidefont = font(14),
            tickfont  = font(11), legendfont = font(12))
end

function downsample(v, n_max)
    n = length(v)
    n <= n_max && return v, 1:n
    step = ceil(Int, n / n_max)
    idx  = 1:step:n
    return v[idx], idx
end

function fit_c(t_cyc, lam)
    n  = length(t_cyc)
    lo = n ÷ 4
    hi = 3 * n ÷ 4
    mean(lam[lo:hi] .* t_cyc[lo:hi] ./ log.(t_cyc[lo:hi]))
end

function main()
    if isempty(ARGS)
        println("Usage: julia scripts/plot_ftmle.jl <ftmle.jld2> [entropy.jld2]")
        return
    end

    ftmle_path = ARGS[1]
    ent_path   = get(ARGS, 2,
        "results/data/backups/backup_data/periodic/entropy_data_springs_alpha_low_periodic.jld2")

    # ── Cargar ftMLE ────────────────────────────────────────────────────────
    d_f     = load(ftmle_path)
    t_cyc   = Float64.(d_f["t_cycles"])
    lam     = Float64.(d_f["lambda"])
    delta_k = Float64(d_f["delta_k"])

    valid  = findall(x -> x > 0 && isfinite(x), lam)
    t_lam  = t_cyc[valid]
    l_lam  = lam[valid]

    # ── Cargar entropía ─────────────────────────────────────────────────────
    d_e       = load(ent_path)
    ent_entry = nothing
    for r in d_e["results"]
        _, dval, _, _ = r
        isapprox(dval, delta_k; atol=1e-10) && (ent_entry = r; break)
    end
    ent_entry === nothing && error("No entry with Δκ=$delta_k in $ent_path")
    _, _, t_ent_raw, S_ent_raw = ent_entry

    # Limitar entropía al rango temporal del ftMLE, excluir t=0
    t_lo  = max(t_lam[1], 1.0)
    t_hi  = t_lam[end]
    mask  = findall(t -> t >= t_lo && t <= t_hi, t_ent_raw)
    t_ent = t_ent_raw[mask]
    S_ent = S_ent_raw[mask]

    # ── Referencia c·ln(t)/t ────────────────────────────────────────────────
    c_fit = fit_c(t_lam, l_lam)
    t_ref = exp10.(range(log10(t_lo), log10(t_hi); length=400))
    l_ref = c_fit .* log.(t_ref) ./ t_ref

    # ── R(t) = λ(t) / [c·ln(t)/t] evaluado en los mismos puntos que λ ───────
    l_ref_at_lam = c_fit .* log.(t_lam) ./ t_lam
    R             = l_lam ./ l_ref_at_lam

    # ── Downsample ───────────────────────────────────────────────────────────
    t_ent_ds, idx_e = downsample(t_ent, 3000)
    S_ent_ds        = S_ent[idx_e]
    t_lam_ds, idx_l = downsample(t_lam, 3000)
    l_lam_ds        = l_lam[idx_l]
    R_ds            = R[idx_l]

    apply_global_plot_style!()

    # ── Panel 1: S̄(t) ───────────────────────────────────────────────────────
    p1 = plot(t_ent_ds, S_ent_ds;
        xscale     = :log10,
        xlabel     = "",
        ylabel     = L"$\bar{S}(t)$",
        label      = latexstring("\\Delta\\kappa = $delta_k"),
        legend     = :topleft,
        framestyle = :box,
        grid       = false,
        lw         = 2.0,
        color      = :navy,
        title      = latexstring("\\Delta\\kappa=$(delta_k),\\ \\alpha=0.1,\\ N=64,\\ \\mathrm{PBC}"),
    )

    # ── Panel 2: λ(t) log-log ───────────────────────────────────────────────
    p2 = plot(t_lam_ds, l_lam_ds;
        xscale     = :log10,
        yscale     = :log10,
        xlabel     = "",
        ylabel     = L"\lambda(t)",
        label      = L"\lambda(t)",
        legend     = :topright,
        framestyle = :box,
        grid       = false,
        lw         = 2.0,
        color      = :crimson,
    )
    plot!(p2, t_ref, l_ref;
        color     = :gray,
        linestyle = :dash,
        lw        = 1.5,
        label     = L"c \cdot \ln(t)/t",
    )

    # ── Panel 3: R(t) = λ(t) / [c·ln(t)/t], semi-log ───────────────────────
    p3 = plot(t_lam_ds, R_ds;
        xscale     = :log10,
        xlabel     = L"t\ \mathrm{(cycles)}",
        ylabel     = L"R(t)",
        label      = false,
        framestyle = :box,
        grid       = false,
        lw         = 2.0,
        color      = :darkgreen,
    )
    hline!(p3, [1.0];
        color     = :gray,
        linestyle = :dash,
        lw        = 1.0,
        label     = false,
    )

    fig = plot(p1, p2, p3; layout=(3,1), size=(800, 950))

    outbase = replace(ftmle_path, ".jld2" => "_fig")
    savefig(fig, outbase * ".pdf")
    println("PDF: $(outbase).pdf")
    savefig(fig, outbase * ".png")
    println("PNG: $(outbase).png")
end

main()
