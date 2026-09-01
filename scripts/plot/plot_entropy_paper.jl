"""
    plot_entropy_paper.jl

Curated entropy figures for the paper: one figure per boundary condition, fixed
(Δκ, color, linestyle, linewidth) tuples from `SPEC_FBC`/`SPEC_PBC`/
`SPEC_PBC_COMPLETE` instead of the cyclic palette `plot_entropy_param_sweep.jl`
uses for exploratory sweeps.

  FBC          : Δκ = 0.05, 0.1, 0.5, 0.7
  PBC          : Δκ = 0.05, 0.1, 0.2, 0.5, 0.7
  PBC complete : Δκ = 0.05 … 0.9, merging the production sweep with extra runs
                 (only produced when more than one PBC file is given)

Usage:
  julia --project=. scripts/plot/plot_entropy_paper.jl <fbc.jld2> <pbc.jld2> [pbc_extra.jld2 ...]
"""

using LaTeXStrings
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils
include("../../src/plot_style.jl");     using .PlotStyle

const USAGE = "Usage: julia --project=. scripts/plot/plot_entropy_paper.jl <fbc.jld2> <pbc.jld2> [pbc_extra.jld2 ...]"
const OUTDIR = "results/figures/entropy"

"""
    paper_curves(paths, spec) -> Vector{Curve}

S̄(t) for each (Δκ, color, linestyle, linewidth) in `spec`, taking the first
result matching Δκ across `paths` — lets a base sweep be extended with extra
runs without overriding it. Missing Δκ are reported and skipped.
"""
function paper_curves(paths, spec)
    # load_results_by_delta merges all files, keyed by Δκ; the dict preserves the first
    # occurrence, so a base file's data takes precedence over extended extras.
    by_delta = load_results_by_delta(paths...)
    curves = Curve[]
    for (d, color, ls, lw) in spec
        haskey(by_delta, d) || (println("  ⚠ Δκ=$d not found"); continue)
        t, E = prepare_ts(by_delta[d])
        t, S = logdownsample(t, entropy_series(E))
        push!(curves, Curve(t, S, latexstring("\\Delta\\kappa = $d");
                            color = color, linestyle = ls, linewidth = lw))
        println("  ✓ Δκ=$d ($(length(t)) points)")
    end
    curves
end

fig(paths, spec) = draw!(logplot(ylabel = L"\bar{S}(t)"), paper_curves(paths, spec))

function main()
    length(ARGS) >= 2 || error(USAGE)
    fbc_path, pbc_paths = ARGS[1], ARGS[2:end]
    apply_style!()

    println("\nFBC ($(basename(fbc_path))):")
    save_fig(fig([fbc_path], SPEC_FBC), OUTDIR, "fig_entropy_FBC")

    println("\nPBC ($(basename(pbc_paths[1]))):")
    # pbc_paths[1:1] extracts only the first PBC file; prevents accidental merging if extra files exist.
    save_fig(fig(pbc_paths[1:1], SPEC_PBC), OUTDIR, "fig_entropy_PBC")

    if length(pbc_paths) > 1
        println("\nPBC complete (Δκ = 0.05–0.9, $(length(pbc_paths)) files):")
        save_fig(fig(pbc_paths, SPEC_PBC_COMPLETE), OUTDIR, "fig_entropy_PBC_complete")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
