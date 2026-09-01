"""
    plot_entropy_param_sweep.jl

Base entropy figure: S̄(t) on logarithmic time, one panel per sweep parameter,
one curve per Δ.

This is the reference shape for every time-series figure in the repo. A new one
is this script with a different `*_curves` function: the loading, styling,
panel and saving are already generic.

Usage:
  julia --project=. scripts/plot/plot_entropy_param_sweep.jl <plot_config.toml | results.jld2>
"""

using JLD2, LaTeXStrings
include("../../src/config.jl");         using .Config
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils
include("../../src/plot_style.jl");     using .PlotStyle

const USAGE = "Usage: julia --project=. scripts/plot/plot_entropy_param_sweep.jl <plot_config.toml | results.jld2>"

plot_config(arg::AbstractString) =
    endswith(lowercase(arg), ".jld2") ? Config.default_plot_config(arg) :
                                        Config.load_plot_config(arg)

"""Results passing the config's param/Δ filters; an empty filter keeps everything."""
function select(results, pcfg)
    keep(v, allowed) = isempty(allowed) ||
                       any(a -> isapprox(Float64(v), a; rtol = 1e-6), allowed)
    filter(r -> keep(r.param, pcfg.filter_params) &&
                keep(r.Delta, pcfg.filter_deltas), results)
end

"""
    curve_label(cfg, res, show_N) -> LaTeXString

Δκ (or Δm) alone, plus N when the file spans several sizes. `compute_trajectories.jl`
can emit a multi-N sweep in one file, and without N the legend would repeat the
same Δ label once per size.
"""
function curve_label(cfg, res, show_N::Bool)
    sym = cfg.system_type == :springs ? "\\Delta \\kappa" : "\\Delta m"
    show_N && hasproperty(res, :N) ?
        latexstring("N = $(res.N),\\ $sym = $(res.Delta)") :
        latexstring("$sym = $(res.Delta)")
end

"""
    entropy_curves(runs, cfg, smooth_delta) -> Vector{Curve}

S̄(t) for each run, styled cyclically. Runs with no positive timestamp cannot go
on a log axis and are reported and dropped.
"""
function entropy_curves(runs, cfg, smooth_delta)
    show_N = length(unique(hasproperty(r, :N) ? r.N : 0 for r in runs)) > 1
    curves = Curve[]
    for (j, res) in enumerate(runs)
        t, E = prepare_ts(res)
        if isempty(t)
            println("    · Δ=$(res.Delta): no positive times, skipped")
            continue
        end
        t, S = logdownsample(t, entropy_series(E; delta = smooth_delta))
        push!(curves, Curve(t, S, curve_label(cfg, res, show_N); cyclic(j)...))
        println("    ✓ Δ=$(res.Delta)$(show_N ? " N=$(res.N)" : "") ($(length(t)) points)")
    end
    curves
end

function main()
    isempty(ARGS) && error(USAGE)
    pcfg = plot_config(ARGS[1])
    apply_style!()

    println("Loading: $(pcfg.input_file)")
    data    = load(pcfg.input_file)
    cfg     = Config.as_experiment_config(data["config"])
    results = select(data["results"], pcfg)

    for p in sort(unique(Float64(r.param) for r in results))
        println("$(cfg.nonlinear) = $p:")
        runs = sort(filter(r -> Float64(r.param) == p, results),
                    by = r -> (hasproperty(r, :N) ? r.N : 0, Float64(r.Delta)))
        fig  = logplot(ylabel = L"\bar{S}(t)")
        draw!(fig, entropy_curves(runs, cfg, pcfg.smooth_delta))
        save_fig(fig, pcfg.outdir_entropy,
                 "entropy_$(cfg.nonlinear)_p$(p)_$(cfg.boundary)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
