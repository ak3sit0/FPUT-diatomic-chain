"""
    plot_entropy_size_sweep.jl

Entropy figure comparing system sizes: S̄(t) for the same (param, Δκ) taken
from several sweep files, one curve per N. Complements
`plot_entropy_param_sweep.jl` (many Δκ, one N) — here it's many N, one Δκ.

N is read from `result.N` when the file stores it per run (size-sweep files
like `nsweep_rise_N*`); single-N sweep files (e.g. the N=64 production sweep)
don't carry it per result, so it falls back to the sweep's own config TOML.

Usage:
  julia --project=. scripts/plot/plot_entropy_size_sweep.jl --param P --delta D <file1.jld2> [file2.jld2 ...]
"""

using JLD2, LaTeXStrings, TOML
include("../../src/config.jl");         using .Config
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils
include("../../src/plot_style.jl");     using .PlotStyle

const USAGE = "Usage: julia --project=. scripts/plot/plot_entropy_size_sweep.jl --param P --delta D <file1.jld2> [file2.jld2 ...]"

function parse_args(args)
    param = delta = nothing
    paths = String[]
    i = 1
    while i <= length(args)
        if args[i] == "--param"
            param = parse(Float64, args[i+1]); i += 2
        elseif args[i] == "--delta"
            delta = parse(Float64, args[i+1]); i += 2
        else
            push!(paths, args[i]); i += 1
        end
    end
    (isnothing(param) || isnothing(delta) || isempty(paths)) && error(USAGE)
    param, delta, paths
end

"""
    result_N(res, config_path) -> Int

System size for a result: its own `N` field when present, otherwise the `N`
recorded in the sweep's config TOML (single-N sweep files fix N for every run
and don't repeat it per result).
"""
function result_N(res, config_path)
    hasproperty(res, :N) && return res.N
    resolved = Config.resolve_config_path(config_path)
    isnothing(resolved) && error("Cannot determine N: no per-result field and config path unresolved ($config_path)")
    TOML.parsefile(resolved)["physics"]["N"]
end

"""S̄(t) for the (param, Δκ) match in each file, one curve per N."""
function size_curves(paths, param, delta)
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
        t, S = logdownsample(t, entropy_series(E))
        push!(curves, Curve(t, S, latexstring("N = $(m.N)"); cyclic(j)...))
        println("  ✓ N=$(m.N) ($(length(t)) points)")
    end
    curves
end

function main()
    param, delta, paths = parse_args(ARGS)
    apply_style!()

    fig = logplot(ylabel = L"\bar{S}(t)")
    draw!(fig, size_curves(paths, param, delta))
    save_fig(fig, "results/figures/entropy",
             "entropy_size_sweep_p$(param)_d$(delta)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
