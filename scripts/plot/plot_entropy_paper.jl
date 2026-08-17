"""
    plot_entropy_paper.jl

Figuras de entropía S̄(t) para el paper, con codificación de color/estilo/grosor
por Δκ (ver results/figures/entropy/instrucciones_figuras_entropia.md).

  FBC          : Δκ = 0.05, 0.1, 0.5, 0.7
  PBC          : Δκ = 0.05, 0.1, 0.2, 0.5, 0.7
  PBC completo : Δκ = 0.05 … 0.9, fusionando el sweep de producción con
                 corridas adicionales (una figura extra, sólo si se pasan
                 más archivos PBC)

Usage:
  julia --project=. scripts/plot/plot_entropy_paper.jl <fbc.jld2> <pbc.jld2> [pbc_extra.jld2 ...]
"""

using JLD2, Plots, LaTeXStrings, Statistics
import Plots: mm
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils

apply_style!() = default(titlefont=font(16), guidefont=font(14),
                         tickfont=font(12),  legendfont=font(12))

"""
    build_fig(jld2_paths, spec, label) -> Plots.Plot

S̄(t) en log-x para cada Δκ de `spec`, tomando el primer resultado que coincide
entre todos los `jld2_paths` (permite extender un sweep con corridas nuevas).
"""
function build_fig(jld2_paths, spec, label)
    by_delta = load_results_by_delta(jld2_paths...)

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
            println("  ⚠ Δκ=$delta no encontrado en $label")
            continue
        end
        res  = by_delta[first(matching)]
        t, E = prepare_ts(res)
        t_ds, S_ds = logdownsample(t, entropy_series(E))
        plot!(fig, t_ds, S_ds;
            color=color, linestyle=ls, lw=lw, alpha=0.9,
            label=latexstring("\\Delta\\kappa = $delta"))
        println("  ✓ Δκ=$delta  ($(length(t)) puntos)")
    end
    fig
end

function save_fig(fig, outdir, stem)
    for ext in (".pdf", ".png")
        savefig(fig, joinpath(outdir, stem * ext))
    end
    println("  → $stem.pdf/.png")
end

function main()
    length(ARGS) >= 2 || error(
        "Usage: julia plot_entropy_paper.jl <fbc.jld2> <pbc.jld2> [pbc_extra.jld2 ...]")
    fbc_path  = ARGS[1]
    pbc_paths = ARGS[2:end]

    apply_style!()
    outdir = "results/figures/entropy"
    mkpath(outdir)

    println("\nFBC (", basename(fbc_path), "):")
    save_fig(build_fig([fbc_path], SPEC_FBC, "FBC"), outdir, "fig_entropy_FBC")

    println("\nPBC (", basename(pbc_paths[1]), "):")
    save_fig(build_fig(pbc_paths[1:1], SPEC_PBC, "PBC"), outdir, "fig_entropy_PBC")

    # Figura extendida sólo cuando hay corridas PBC adicionales que fusionar.
    if length(pbc_paths) > 1
        println("\nPBC completo (Δκ = 0.05–0.9, ", length(pbc_paths), " archivos):")
        save_fig(build_fig(pbc_paths, SPEC_PBC_COMPLETE, "PBC-complete"),
                 outdir, "fig_entropy_PBC_complete")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
