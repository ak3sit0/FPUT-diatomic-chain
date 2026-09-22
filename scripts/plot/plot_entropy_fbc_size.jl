"""
    plot_entropy_fbc_size.jl

The paper's entropy figure (Δκ curves in the `SPEC_FBC` / `SPEC_PBC` styles) for
any system size and either boundary condition. `plot_entropy_paper.jl` builds
the same figures but for N = 64 only (it takes one file per boundary and no size
selector); this one merges several sweep files, keeps the results of one N, and
names the output after it, so the N = 128 and N = 256 versions sit next to the
N = 64 original.

  --bc fbc   Δκ = 0.05, 0.1, 0.5, 0.7           (first-mode excitation)
  --bc pbc   Δκ = 0.05, 0.1, 0.2, 0.5, 0.7      (mode-2 excitation)
  --N n      keep only results of this size (needed when a file mixes sizes)
  --tmin t   left edge of the time axis (default: the first decade with data)
  --norm     raw (default) or logN

S̄(t) is the raw time-averaged spectral entropy, as in the paper figure, so its
ceiling is ln N and moves with N; `--norm logN` divides by it for a cross-N
comparison. The entropy is computed on the full saved series before the axis is
cut, since the growing-window average depends on the whole history.

When a Δκ appears in several files the first one wins, so list the preferred
file first.

Usage:
  julia --project=. scripts/plot/plot_entropy_fbc_size.jl [--bc fbc|pbc] [--N n] [--tmin t] [--norm raw|logN] <a.jld2> [b.jld2 ...]
"""

using JLD2, LaTeXStrings, TOML
include("../../src/config.jl");         using .Config
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils
include("../../src/plot_style.jl");     using .PlotStyle

const USAGE = "Usage: julia --project=. scripts/plot/plot_entropy_fbc_size.jl [--bc fbc|pbc] [--N n] [--tmin t] [--norm raw|logN] <a.jld2> [b.jld2 ...]"
const OUTDIR = "results/figures/entropy/paper"
const SPECS = Dict("fbc" => SPEC_FBC, "pbc" => SPEC_PBC)

function parse_args(args)
    o = (bc = "fbc", N = nothing, tmin = nothing, norm = "raw")
    paths, i = String[], 1
    while i <= length(args)
        a = args[i]
        if     a == "--bc";   o = merge(o, (bc = args[i+1],));                    i += 2
        elseif a == "--N";    o = merge(o, (N = parse(Int, args[i+1]),));         i += 2
        elseif a == "--tmin"; o = merge(o, (tmin = parse(Float64, args[i+1]),));  i += 2
        elseif a == "--norm"; o = merge(o, (norm = args[i+1],));                  i += 2
        else push!(paths, a); i += 1
        end
    end
    (isempty(paths) || !haskey(SPECS, o.bc) || o.norm ∉ ("raw", "logN")) && error(USAGE)
    o, paths
end

"""System size of a sweep file: the per-result `N` if stored, else the config TOML."""
function file_N(d)
    res = first(d["results"])
    hasproperty(res, :N) && return res.N
    resolved = Config.resolve_config_path(d["config"])
    isnothing(resolved) && error("Cannot determine N: results carry no N and the config is unresolved")
    TOML.parsefile(resolved)["physics"]["N"]
end

"""Results of size `N` across `paths`, keyed by Δκ, first occurrence wins."""
function results_by_delta(paths, N)
    by_delta = Dict{Float64,Any}()
    for p in paths
        d = load(p)
        n_file = file_N(d)
        for r in d["results"]
            n = hasproperty(r, :N) ? r.N : n_file
            n == N || continue
            haskey(by_delta, Float64(r.Delta)) || (by_delta[Float64(r.Delta)] = r)
        end
    end
    by_delta
end

function main()
    o, paths = parse_args(ARGS)
    apply_style!()
    N = isnothing(o.N) ? file_N(load(first(paths))) : o.N
    scale = o.norm == "logN" ? log(N) : 1.0

    by_delta = results_by_delta(paths, N)
    isempty(by_delta) && error("No results with N = $N in: $(join(paths, ", "))")

    curves = Curve[]
    for (d, color, ls, lw) in SPECS[o.bc]
        haskey(by_delta, d) || (println("  ⚠ Δκ=$d not found for N=$N"); continue)
        t, E = prepare_ts(by_delta[d])
        t, S = logdownsample(t, entropy_series(E) ./ scale)
        push!(curves, Curve(t, S, latexstring("\\Delta\\kappa = $d");
                            color = color, linestyle = ls, linewidth = lw))
        println("  ✓ N=$N Δκ=$d  ($(length(t)) points, t=[$(round(t[1]; sigdigits=3)), $(round(t[end]; sigdigits=3))], S_end=$(round(S[end]; digits=3)))")
    end
    isempty(curves) && error("None of the $(uppercase(o.bc)) Δκ values found for N = $N")

    tmax = maximum(last(c.x) for c in curves)
    lo = isnothing(o.tmin) ? exp10(ceil(log10(maximum(first(c.x) for c in curves)))) : o.tmin
    hi = exp10(ceil(log10(tmax)))
    ylabel = o.norm == "logN" ? L"\bar{S}(t)\,/\,\ln N" : L"\bar{S}(t)"
    fig = logplot(; ylabel = ylabel, xlims = (lo, hi),
                  decades = Int(round(log10(lo))):Int(round(log10(hi))))
    draw!(fig, curves)
    stem = "fig_entropy_$(uppercase(o.bc))_N$(N)" * (o.norm == "logN" ? "_logN" : "")
    save_fig(fig, OUTDIR, stem)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
