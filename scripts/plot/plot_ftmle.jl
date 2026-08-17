"""
    plot_ftmle.jl

Figura de un panel: λ(t) en escala log-log + ley de potencias ajustada t^{-δ}.
El ajuste se hace sobre el 50% central de la curva del caso Δκ=0.1 (representativo).

Soporta uno o varios archivos JLD2 de ftMLE superpuestos.
Cada caso recibe un color y estilo de línea distintos.

Usage:
  julia --project=. scripts/plot_ftmle.jl <ftmle1.jld2> [ftmle2.jld2 ...]
"""

using JLD2, Plots, LaTeXStrings, Statistics, Dates
import Plots: mm
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils

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

"""
    fit_powerlaw(t, lam; t_max, n_bins) -> (c, δ)

Binning logarítmico + regresión lineal en log-log.
Divide [t[1], t_max] en n_bins bins de igual ancho en log(t),
promedia λ dentro de cada bin, y ajusta sobre las medias de bin.
Cada década pesa igual independientemente de la densidad de puntos.
"""
function fit_powerlaw(t, lam; t_max=1e4, n_bins=20)
    mask = (t .<= t_max) .& (t .> 0) .& (lam .> 0)
    t_m  = t[mask]
    l_m  = lam[mask]
    isempty(t_m) && error("Sin datos en el rango de ajuste")

    edges = exp10.(range(log10(t_m[1]), log10(t_max); length=n_bins+1))
    t_bin = Float64[]; l_bin = Float64[]
    for i in 1:n_bins
        idx = findall(s -> edges[i] <= s < edges[i+1], t_m)
        isempty(idx) && continue
        push!(t_bin, exp(mean(log.(t_m[idx]))))   # media geométrica del bin
        push!(l_bin, exp(mean(log.(l_m[idx]))))   # media geométrica de λ
    end

    x  = log.(t_bin); y = log.(l_bin)
    mx = mean(x);     my = mean(y)
    slope     = sum((x .- mx) .* (y .- my)) / sum((x .- mx).^2)
    intercept = my - slope * mx
    return exp(intercept), -slope
end

const DELTA_REF = 0.05   # caso del que se toma la referencia

function main()
    if isempty(ARGS)
        println("Usage: julia scripts/plot_ftmle.jl <ftmle1.jld2> [ftmle2.jld2 ...]")
        return
    end

    ftmle_paths = ARGS[:]
    apply_global_plot_style!()

    fig = plot(; xscale=:log10, yscale=:log10,
                 ylabel     = L"\lambda(t)",
                 xlabel     = L"t\ \mathrm{(cycles)}",
                 legend     = :topright,
                 framestyle = :box,
                 grid       = true,
                 gridalpha  = 0.2,
                 gridstyle  = :dot,
                 left_margin   = 5mm,
                 right_margin  = 4mm,
                 top_margin    = 4mm,
                 bottom_margin = 6mm,
                 size          = (800, 500))

    ref_c  = nothing
    ref_δ  = nothing
    ref_t0 = nothing
    ref_t1 = nothing

    for (k, ftmle_path) in enumerate(ftmle_paths)
        color = cyc(FTMLE_PALETTE, k)
        ls    = cyc(FTMLE_LINESTYLES, k)

        d_f     = load(ftmle_path)
        t_cyc   = Float64.(d_f["t_cycles"])
        lam     = Float64.(d_f["lambda"])
        delta_k = Float64(d_f["delta_k"])

        valid = findall(x -> x > 0 && isfinite(x), lam)
        t_lam = t_cyc[valid]
        l_lam = lam[valid]

        if isapprox(delta_k, DELTA_REF; atol=1e-10)
            ref_c, ref_δ = fit_powerlaw(t_lam, l_lam)
            ref_t0 = t_lam[1]
            ref_t1 = t_lam[end]
        end

        lw = isapprox(delta_k, DELTA_REF; atol=1e-10) ? 3.5 : 1.5
        t_lam_ds, idx_l = downsample(t_lam, 3000)
        plot!(fig, t_lam_ds, l_lam[idx_l];
            lw=lw, color=color, linestyle=ls,
            label=latexstring("\\Delta\\kappa = $(delta_k)"))
    end

    # ── Referencia c·t^{-δ} sobre la segunda mitad ──────────────────────────
    if !isnothing(ref_c)
        t_ref   = exp10.(range(log10(ref_t0), log10(ref_t1); length=400))
        δ_str   = string(round(ref_δ; digits=3))
        plot!(fig, t_ref, ref_c .* t_ref .^ (-ref_δ);
            color=:black, linestyle=:solid, lw=2.0,
            label=latexstring("c\\,t^{-$(δ_str)}"))
        println("Ajuste (Δκ=$(DELTA_REF)):  c = $(round(ref_c; sigdigits=3))   δ = $(round(ref_δ; digits=4))")
    end

    figdir = "results/figures/liapunov_exponent"
    mkpath(figdir)
    # Detectar boundary: usar el campo del JLD2 si existe, si no asumir "periodic"
    bc_tags = unique([get(load(p), "boundary", "periodic") for p in ftmle_paths])
    bc_str  = length(bc_tags) == 1 ? bc_tags[1] : "mixed"
    tag     = join(replace.(string.(sort(unique(
                  [Float64(load(p)["delta_k"]) for p in ftmle_paths]))), "." => "p"), "_")
    outbase = joinpath(figdir, "ftmle_$(bc_str)_delta$(tag)_$(today())")
    savefig(fig, outbase * ".pdf")
    println("PDF: $(outbase).pdf")
end

main()
