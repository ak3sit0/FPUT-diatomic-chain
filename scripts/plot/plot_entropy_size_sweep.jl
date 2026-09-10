"""
    plot_entropy_size_sweep.jl

Entropy figure comparing system sizes: S̄(t) for the same (param, Δκ) taken
from several sweep files, one curve per N. Complements
`plot_entropy_param_sweep.jl` (many Δκ, one N) — here it's many N, one Δκ.

Raw S̄ is not comparable across N: it lives in [0, ln N], so the same curve
height means a different degree of spreading at N=64 and N=256. Every curve is
therefore normalized before plotting (`--norm`):

  logN  S̄/ln N — fraction of the equipartition ceiling, 1 ⇔ every mode carries
        the same energy. Bounded in [0,1] and directly comparable.
  neff  e^{S̄} — effective number of excited modes. Use this to tell a packet
        that spreads over a *fixed number of modes* (n_eff flat in N) from one
        that spreads over a *fixed fraction of the Brillouin zone* (n_eff ∝ N):
        S̄/ln N shrinks with N in the first case purely through its denominator.
  raw   S̄ — escape hatch; only meaningful at fixed N.

N is read from `result.N` when the file stores it per run (size-sweep files
like `nsweep_rise_N*`); single-N sweep files (e.g. the N=64 production sweep)
don't carry it per result, so it falls back to the sweep's own config TOML.

Usage:
  julia --project=. scripts/plot/plot_entropy_size_sweep.jl --param P --delta D [--norm logN|neff|raw] <file1.jld2> [file2.jld2 ...]
"""

using JLD2, LaTeXStrings, TOML, Printf
include("../../src/config.jl");         using .Config
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils
include("../../src/plot_style.jl");     using .PlotStyle

const USAGE = "Usage: julia --project=. scripts/plot/plot_entropy_size_sweep.jl --param P --delta D [--norm logN|neff|raw] <file1.jld2> [file2.jld2 ...]"

"""
Normalization of S̄ for cross-N comparison: how to map `(S, N)` to the plotted
ordinate, how to label it, the axis settings it implies, and the reference level
worth drawing (`nothing` when the quantity has no N-independent ceiling).
"""
struct Normalization
    apply::Function
    ylabel::LaTeXString
    axis::NamedTuple
    ceiling::Union{Float64,Nothing}
end

const NORMALIZATIONS = Dict(
    # Bounded in [0,1]: pin the axis there so figures at different Δκ are
    # visually comparable and the equipartition ceiling keeps its meaning.
    "logN" => Normalization((S, N) -> S ./ log(N),
                            L"\bar{S}(t)\,/\,\ln N",
                            (ylims = (0.0, 1.02), yticks = 0:0.2:1.0),
                            1.0),
    # Spans 1 → N over the run, so log-y keeps the early spreading visible.
    "neff" => Normalization((S, N) -> exp.(S),
                            L"n_{\mathrm{eff}}(t) = \exp[\,\bar{S}(t)\,]",
                            (yscale = :log10,),
                            nothing),
    "raw"  => Normalization((S, N) -> S,
                            L"\bar{S}(t)",
                            NamedTuple(),
                            nothing),
)

function parse_args(args)
    param = delta = nothing
    norm = "logN"
    paths = String[]
    i = 1
    while i <= length(args)
        if args[i] == "--param"
            param = parse(Float64, args[i+1]); i += 2  # skip both the flag and its value
        elseif args[i] == "--delta"
            delta = parse(Float64, args[i+1]); i += 2
        elseif args[i] == "--norm"
            norm = args[i+1]; i += 2
        else
            push!(paths, args[i]); i += 1  # positional argument ⇒ a file path
        end
    end
    (isnothing(param) || isnothing(delta) || isempty(paths)) && error(USAGE)
    haskey(NORMALIZATIONS, norm) ||
        error("Unknown --norm '$norm'. Valid: " * join(sort(collect(keys(NORMALIZATIONS))), ", "))
    param, delta, norm, paths
end

"""
    result_N(res, config_path) -> Int

System size for a result: its own `N` field when present, otherwise the `N`
recorded in the sweep's config TOML (single-N sweep files fix N for every run
and don't repeat it per result).
"""
function result_N(res, config_path)
    hasproperty(res, :N) && return res.N  # per-result N (multi-N sweeps) takes precedence
    # Fallback: parse the config TOML for a single fixed N recorded in every sweep file.
    # This branch is taken only when results have no :N field (single-N production sweeps).
    resolved = Config.resolve_config_path(config_path)
    isnothing(resolved) && error("Cannot determine N: no per-result field and config path unresolved ($config_path)")
    TOML.parsefile(resolved)["physics"]["N"]
end

"""Normalized S̄(t) for the (param, Δκ) match in each file, one curve per N."""
function size_curves(paths, param, delta, norm::Normalization)
    matches = NamedTuple[]
    for path in paths
        d = load(path)
        idx = findfirst(r -> isapprox(Float64(r.param), param; rtol=1e-6) &&
                              isapprox(Float64(r.Delta), delta; rtol=1e-6), d["results"])
        if isnothing(idx)
            println("  ⚠ param=$param, Δκ=$delta not found in $(basename(path))")
            continue
        end
        res = d["results"][idx]
        push!(matches, (N = result_N(res, d["config"]), res = res))
    end
    sort!(matches, by = m -> m.N)

    curves = Curve[]
    for (j, m) in enumerate(matches)
        t, E = prepare_ts(m.res)
        # Entropy on the full series, then downsample: the growing-window average
        # is a fixed fraction of elapsed time, so thinning first would change it.
        S = entropy_series(E)
        t, y = logdownsample(t, norm.apply(S, m.N))
        push!(curves, Curve(t, y, latexstring("N = $(m.N)"); cyclic(j; lw = 2.5)...))
        @printf("  ✓ N=%-4d %5d pts   t_max=%.2e   S=%.3f   S/lnN=%.3f   n_eff=%.1f\n",
                m.N, length(t), t[end], S[end], S[end]/log(m.N), exp(S[end]))
    end
    curves
end

function main()
    param, delta, norm_key, paths = parse_args(ARGS)
    norm = NORMALIZATIONS[norm_key]
    apply_style!()

    curves = size_curves(paths, param, delta, norm)
    isempty(curves) && error("No curve matched param=$param, Δκ=$delta in any input file")

    # Start the axis where *every* curve has data. Runs saved at different
    # cadences (the N=64 sweep stores 200k samples, the size sweeps 4k), so the
    # latest-starting curve sets the floor — ceil (not floor) of its first point,
    # otherwise the opening decade still shows only the earlier curves alone.
    lo = exp10(ceil(log10(maximum(first(c.x) for c in curves))))
    hi = exp10(ceil(log10(maximum(last(c.x) for c in curves))))
    fig = logplot(; ylabel = norm.ylabel, xlims = (lo, hi),
                  decades = Int(log10(lo)):Int(log10(hi)), norm.axis...)
    draw!(fig, curves)
    isnothing(norm.ceiling) ||
        guide_hline!(fig, norm.ceiling; label = "equipartition")
    save_fig(fig, "results/figures/entropy/size_sweep",
             "entropy_size_sweep_p$(param)_d$(delta)_$(norm_key)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
