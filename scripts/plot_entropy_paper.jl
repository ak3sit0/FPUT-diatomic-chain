"""
    plot_entropy_paper.jl

Figuras de entropía S̄(t) para el paper con codificación de color/grosor
según instrucciones_figuras_entropia.md.

FBC: 4 curvas (Δκ = 0.05, 0.1, 0.5, 0.7)
PBC: 5 curvas (Δκ = 0.05, 0.1, 0.2, 0.5, 0.7)

Usage:
  julia --project=. scripts/plot_entropy_paper.jl <fbc.jld2> <pbc.jld2>
"""

using JLD2, Plots, LaTeXStrings, Statistics
import Plots: mm
include("../src/fput_analysis.jl"); using .FPUTAnalysis

const SMOOTH_DELTA = 0.6
const EPS          = 1e-18

# ── Spec: (delta, color_hex, linestyle, lw) ──────────────────────────────────
const SPEC_FBC = [
    (0.05, "#1a3a6b", :solid,   4.0),
    (0.1,  "#2171b5", :dash,    2.5),
    (0.5,  "#d62728", :dashdot, 2.0),
    (0.7,  "#f4845f", :dot,     1.5),
]

const SPEC_PBC = [
    (0.05, "#1a3a6b", :solid,   4.0),
    (0.1,  "#2171b5", :dash,    2.5),
    (0.2,  "#6baed6", :dot,     2.5),
    (0.5,  "#d62728", :dashdot, 2.0),
    (0.7,  "#f4845f", :solid,   1.5),
]

function apply_style!()
    default(titlefont=font(16), guidefont=font(14),
            tickfont=font(12),  legendfont=font(12))
end

function compute_entropy(modal_E::Matrix)
    smoothed = similar(modal_E, Float64)
    for i in 1:size(modal_E, 1)
        smoothed[i, :] = FPUTAnalysis.sliding_window_avg(
            Float64.(modal_E[i, :]), SMOOTH_DELTA)
    end
    total  = sum(smoothed, dims=1)
    p      = smoothed ./ (total .+ EPS)
    p_safe = p .+ (p .== 0.0)
    -vec(sum(p .* log.(p_safe), dims=1))
end

function load_results_by_delta(jld2_path)
    data    = load(jld2_path)
    results = data["results"]
    by_delta = Dict{Float64, Any}()
    for res in results
        d = Float64(res.Delta)
        haskey(by_delta, d) || (by_delta[d] = res)
    end
    by_delta
end

function prepare_ts(res)
    t = Float64.(Vector(res.scaled_t))
    E = Float64.(Matrix(res.modal_E))
    size(E, 1) > size(E, 2) && (E = E')
    idx = findall(>(1e-3), t)
    t[idx], E[:, idx]
end

function logdownsample(t, y, n_max=3000)
    n = length(t)
    n <= n_max && return t, y
    log_edges = range(log10(t[1]), log10(t[end]); length=n_max+1)
    idx = Int[]
    for i in 1:n_max
        lo, hi = 10^log_edges[i], 10^log_edges[i+1]
        j = findfirst(x -> lo <= x < hi, t)
        isnothing(j) || push!(idx, j)
    end
    isempty(idx) && return t, y
    t[idx], y[idx]
end

function build_fig(jld2_path, spec, bc_label)
    by_delta = load_results_by_delta(jld2_path)

    exp_range = 0:6
    fig = plot(;
        xlabel     = L"t",
        ylabel     = L"\bar{S}(t)",
        xscale     = :log10,
        legend     = :bottomright,
        framestyle = :box,
        grid       = true, gridalpha = 0.2, gridstyle = :dot,
        xticks     = (10.0 .^ exp_range,
                      [latexstring("10^{$(i)}") for i in exp_range]),
        left_margin   = 5mm,  right_margin  = 3mm,
        top_margin    = 2mm,  bottom_margin = 5mm,
        size          = (700, 500))

    for (delta, color, ls, lw) in spec
        matching = filter(k -> isapprox(k, delta; atol=1e-9), collect(keys(by_delta)))
        if isempty(matching)
            println("  ⚠ Δκ=$delta no encontrado en $bc_label")
            continue
        end
        res = by_delta[matching[1]]
        t, E = prepare_ts(res)
        S    = compute_entropy(E)
        t_ds, S_ds = logdownsample(t, S)
        plot!(fig, t_ds, S_ds;
            color=color, linestyle=ls, lw=lw, alpha=0.9,
            label=latexstring("\\Delta\\kappa = $delta"))
        println("  ✓ Δκ=$delta  ($(length(t)) puntos)")
    end
    fig
end

function main()
    length(ARGS) == 2 || error(
        "Usage: julia plot_entropy_paper.jl <fbc.jld2> <pbc.jld2>")
    fbc_path = ARGS[1]
    pbc_path = ARGS[2]
    apply_style!()

    outdir = "results/figures/entropy"
    mkpath(outdir)

    println("\nFBC (", basename(fbc_path), "):")
    fig_fbc = build_fig(fbc_path, SPEC_FBC, "FBC")
    savefig(fig_fbc, joinpath(outdir, "fig_entropy_FBC.pdf"))
    savefig(fig_fbc, joinpath(outdir, "fig_entropy_FBC.png"))
    println("  → fig_entropy_FBC.pdf/.png")

    println("\nPBC (", basename(pbc_path), "):")
    fig_pbc = build_fig(pbc_path, SPEC_PBC, "PBC")
    savefig(fig_pbc, joinpath(outdir, "fig_entropy_PBC.pdf"))
    savefig(fig_pbc, joinpath(outdir, "fig_entropy_PBC.png"))
    println("  → fig_entropy_PBC.pdf/.png")
end

main()
