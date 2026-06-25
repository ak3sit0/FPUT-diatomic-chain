"""
    plot_entropy_param_sweep.jl

Plot entropy and localization curves for parameter sweeps using the same
grouping and visual style as the recovery scripts, but adapted to this repo's
JLD2 layout.

Usage:
  julia --project=. scripts/plot_entropy_param_sweep.jl <plot_config.toml | results.jld2>
"""

using Plots, JLD2, LaTeXStrings
include("../src/config.jl"); using Main.Config
include("../src/fput_analysis.jl"); using .FPUTAnalysis

const SMOOTH_DELTA = 0.6
const EPS = 1e-18
const LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

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

function smooth_modal_energies(modal_E::AbstractMatrix, delta::Float64)
    smoothed = similar(modal_E, Float64)
    for i in 1:size(modal_E, 1)
        smoothed[i, :] = FPUTAnalysis.sliding_window_avg(Float64.(modal_E[i, :]), delta)
    end
    smoothed
end

function entropy_from_smooth(E_smooth)
    total_E = sum(E_smooth, dims=1)
    p = E_smooth ./ (total_E .+ EPS)
    p_safe = p .+ (p .== 0.0)
    -vec(sum(p .* log.(p_safe), dims=1))
end

function compute_entropy_smooth(modal_E::Matrix; delta::Float64=SMOOTH_DELTA)
    entropy_from_smooth(smooth_modal_energies(modal_E, delta))
end


function load_modal_results(path::String)
    isfile(path) || error("File not found: $path")
    println("Loading: $path")
    data = load(path)
    results = data["results"]
    config = data["config"]
    if isa(config, String) && isfile(config)
        config = Config.load_experiment_config(config)
    end
    return results, config
end

function prepare_run(res)
    t_full = Float64.(Vector(res.scaled_t))
    E_full = Float64.(Matrix(res.modal_E))
    if size(E_full, 1) > size(E_full, 2)
        E_full = E_full'
    end
    idx = findall(>(1e-3), t_full)
    ProcessedRun(Float64(res.Delta), t_full[idx], E_full[:, idx])
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
    default(titlefont = font(16), guidefont = font(14), tickfont = font(11), legendfont = font(12))
end

function build_base_plot(ylabel::LaTeXString)
    plot(xlabel = L"t", ylabel = ylabel, legend = :bottomright, framestyle = :box,
         grid = false, xscale = :log10, size = (800, 600))
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

        S = compute_entropy_smooth(run.modal_E; delta = delta_smooth)
        style = LINESTYLES[mod1(j, length(LINESTYLES))]
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
