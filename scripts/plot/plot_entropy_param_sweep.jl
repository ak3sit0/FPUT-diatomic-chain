"""
    plot_entropy_param_sweep.jl

Plot entropy and localization curves for parameter sweeps using the same
grouping and visual style as the recovery scripts, but adapted to this repo's
JLD2 layout.

Usage:
  julia --project=. scripts/plot_entropy_param_sweep.jl <plot_config.toml | results.jld2>
"""

using Plots, JLD2, LaTeXStrings
import Plots: mm
include("../../src/config.jl");         using .Config
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils

struct ProcessedRun
    delta::Float64
    times::Vector{Float64}
    modal_E::Matrix{Float64}
end

function build_plot_config(input_arg::Union{String,Nothing})
    if isnothing(input_arg)
        error("Usage: julia scripts/plot_entropy_param_sweep.jl <plot_config.toml | results.jld2>")
    elseif endswith(lowercase(input_arg), ".jld2")
        return Config.default_plot_config(input_arg)
    else
        return Config.load_plot_config(input_arg)
    end
end

function load_modal_results(path::String)
    isfile(path) || error("File not found: $path")
    println("Loading: $path")
    data = load(path)
    return data["results"], Config.as_experiment_config(data["config"])
end

function prepare_run(res)
    t, E = prepare_ts(res)
    ProcessedRun(Float64(res.Delta), t, E)
end

function group_results_by_param(results)
    grouped = Dict{Float64, Vector{ProcessedRun}}()
    for res in results
        param_val = Float64(res.param)
        push!(get!(grouped, param_val, Vector{ProcessedRun}()), prepare_run(res))
    end
    for runs in values(grouped)
        sort!(runs, by = r -> r.delta)
    end
    grouped
end

function filter_results(results, filter_params::Vector{Float64}, filter_deltas::Vector{Float64})
    isempty(filter_params) && isempty(filter_deltas) && return results
    filter(r ->
        (isempty(filter_params) || any(p -> isapprox(Float64(r.param), p, rtol=1e-6), filter_params)) &&
        (isempty(filter_deltas) || any(d -> isapprox(Float64(r.Delta), d, rtol=1e-6), filter_deltas)),
        results)
end

function apply_global_plot_style!()
    default(titlefont = font(16), guidefont = font(18), tickfont = font(13), legendfont = font(13))
end

function build_base_plot(ylabel::LaTeXString)
    exp_range = 0:6
    plot(xlabel = L"t", ylabel = ylabel, legend = :bottomright, framestyle = :box,
         grid = false, xscale = :log10, size = (800, 600),
         xticks = (10.0 .^ exp_range, [latexstring("10^{$(i)}") for i in exp_range]),
         left_margin = 6mm, right_margin = 4mm, top_margin = 2mm, bottom_margin = 6mm)
end

function label_for_delta(config, delta)
    if config.system_type == :springs
        latexstring("\\Delta \\kappa = $delta")
    else
        latexstring("\\Delta m = $delta")
    end
end

function plot_entropy_for_param(param_val, runs, config, delta_smooth, outdir_entropy)
    println("  Parameter $(config.nonlinear)=$param_val:")
    p_S = build_base_plot(L"$\bar{S}(t)$")

    for (j, run) in enumerate(runs)
        if isempty(run.times)
            println("    · Skipping Δ=$(run.delta) (no positive timestamps)")
            continue
        end

        S = entropy_series(run.modal_E; delta = delta_smooth)
        style = cyc(LINESTYLES, j)
        label_str = label_for_delta(config, run.delta)

        plot!(p_S, run.times, S, label = label_str, lw = 2.0, linestyle = style, alpha = 0.8)
        println("    ✓ Δ=$(run.delta) ($(length(run.times)) points)")
    end

    fname_S = joinpath(outdir_entropy, "entropy_$(config.nonlinear)_p$(param_val)_$(config.boundary).pdf")
    savefig(p_S, fname_S)
    println("Saved: $fname_S\n")
end

function main()
    if isempty(ARGS)
        error("Usage: julia scripts/plot_entropy_param_sweep.jl <plot_config.toml | results.jld2>")
    end

    pcfg = build_plot_config(ARGS[1])
    apply_global_plot_style!()
    mkpath(pcfg.outdir_entropy)

    results, config = load_modal_results(pcfg.input_file)
    results = filter_results(results, pcfg.filter_params, pcfg.filter_deltas)
    grouped = group_results_by_param(results)

    println("Processing entropy with moving average...\n")

    for param_val in sort(collect(keys(grouped)))
        plot_entropy_for_param(param_val, grouped[param_val], config,
            pcfg.smooth_delta, pcfg.outdir_entropy)
    end

    println("✓ Done. Entropy plots saved to:")
    println("  • $(pcfg.outdir_entropy)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
