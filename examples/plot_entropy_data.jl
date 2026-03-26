# Load precomputed MODAL ENERGIES and plot entropy with moving average
# Usage: julia --project=. examples/plot_entropy_data.jl configs/plot_entropy_alpha_periodic.toml

using Plots, LaTeXStrings, JLD2, Statistics
include("../src/config.jl");          using Main.Config
include("../src/parameters.jl");      using Main.Parameters
include("../src/energy_analysis.jl"); using Main.EnergyAnalysis

# Style constants (not experiment-specific — change here or in build_plot config)
const SMOOTH_DELTA = 0.6
const EPS = 1e-18
const LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

struct ProcessedRun
    delta::Float64
    times::Vector{Float64}
    modal_E::Matrix{Float64}
end

# ==========================================================================
# PHYSICAL HELPERS
# ==========================================================================

function smooth_modal_energies(modal_E::AbstractMatrix, delta::Float64)
    smoothed = similar(modal_E, Float64)
    #for i in 1:size(modal_E, 1)
    for i in eachindex(modal_E, 1)
        smoothed[i, :] = EnergyAnalysis.sliding_delta_average(Float64.(modal_E[i, :]), delta)
    end
    return smoothed
end

function entropy_from_smooth(E_smooth)
    total_E = sum(E_smooth, dims=1)
    p = E_smooth ./ (total_E .+ EPS)
    p_safe = p .+ (p .== 0.0)
    -vec(sum(p .* log.(p_safe), dims=1))
end

function compute_entropy_smooth(modal_E::Matrix; delta::Float64=SMOOTH_DELTA)
    E_smooth = smooth_modal_energies(modal_E, delta)
    entropy_from_smooth(E_smooth)
end

function compute_xi_from_energies(modal_E::Matrix; delta::Float64=SMOOTH_DELTA)
    E_smooth = smooth_modal_energies(modal_E, delta)
    N = size(E_smooth, 1)
    sum_all = vec(sum(E_smooth, dims=1))
    mid = div(N, 2)

    if mid >= N
        sum_upper = zeros(length(sum_all))
        w_k = zeros(0, length(sum_all))
        S_upper = zeros(length(sum_all))
    else
        E_upper = E_smooth[mid+1:end, :]
        sum_upper = vec(sum(E_upper, dims=1))
        w_k = E_upper ./ (sum_upper .+ EPS)'
        w_k .= max.(w_k, 1e-12)
        S_upper = -vec(sum(w_k .* log.(w_k), dims=1))
    end

    @. (2 * sum_upper / (sum_all .+ EPS)) * exp(S_upper) / N
end

# ==========================================================================
# DATA PREPARATION
# ==========================================================================

function load_modal_results(path::String)
    if !isfile(path)
        error("File not found: $path")
    end
    println("Loading: $path")
    data = load(path)
    (data["results"], data["config"])
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
        sort!(runs, by=r->r.delta)
    end
    grouped
end

"""
    filter_results(results, filter_params, filter_deltas) -> filtered results

Keep only runs whose param ∈ filter_params AND delta ∈ filter_deltas.
If either filter vector is empty, that dimension is not filtered (include all).
"""
function filter_results(results, filter_params::Vector{Float64}, filter_deltas::Vector{Float64})
    isempty(filter_params) && isempty(filter_deltas) && return results
    filter(r ->
        (isempty(filter_params) || any(p -> isapprox(Float64(r.param), p, rtol=1e-6), filter_params)) &&
        (isempty(filter_deltas) || any(d -> isapprox(Float64(r.Delta), d, rtol=1e-6), filter_deltas)),
        results
    )
end

# ==========================================================================
# PLOTTING HELPERS
# ==========================================================================

function apply_global_plot_style!()
    default(titlefont=font(16), guidefont=font(14), tickfont=font(11), legendfont=font(12))
end

function build_base_plot(ylabel::LaTeXString)
    plot(xlabel=L"t", ylabel=ylabel, legend=:bottomright, framestyle=:box,
         grid=false, xscale=:log10, size=(800, 600))
end

function label_for_delta(config, delta)
    if config.system_type == :springs
        latexstring("\\Delta \\kappa = $delta")
    else
        latexstring("\\Delta m = $delta")
    end
end

function plot_entropy_xi_for_param(param_val, runs, config, delta_smooth, outdir_entropy, outdir_xi)
    println("  Parameter $(config.nonlinear)=$param_val:")
    p_S = build_base_plot(L"$\bar{S}(t)$")
    p_xi = build_base_plot(L"$\xi(t)$")

    for (j, run) in enumerate(runs)
        if isempty(run.times)
            println("    · Skipping Δ=$(run.delta) (no positive timestamps)")
            continue
        end

        S = compute_entropy_smooth(run.modal_E; delta=delta_smooth)
        xi = compute_xi_from_energies(run.modal_E; delta=delta_smooth)
        style = LINESTYLES[mod1(j, length(LINESTYLES))]
        label_str = label_for_delta(config, run.delta)

        plot!(p_S, run.times, S, label=label_str, lw=2.0, linestyle=style, alpha=0.8)
        plot!(p_xi, run.times, xi, label=label_str, lw=2.0, linestyle=style, alpha=0.8)
        println("    ✓ Δ=$(run.delta) ($(length(run.times)) points)")
    end

    fname_S = joinpath(outdir_entropy, "entropy_$(config.nonlinear)_p$(param_val)_$(config.boundary).pdf")
    savefig(p_S, fname_S)
    println("Saved: $fname_S")

    fname_xi = joinpath(outdir_xi, "xi_$(config.nonlinear)_p$(param_val)_$(config.boundary).pdf")
    savefig(p_xi, fname_xi)
    println("Saved: $fname_xi\n")
end

# ==========================================================================
# MAIN ENTRY
# ==========================================================================

function main()
    if isempty(ARGS)
        error("Usage: julia plot_entropy_data.jl <plot_config.toml>\n" *
              "Example: julia --project=. examples/plot_entropy_data.jl " *
              "configs/plot_entropy_alpha_periodic.toml")
    end

    pcfg = Config.load_plot_config(ARGS[1])
    apply_global_plot_style!()
    mkpath(pcfg.outdir_entropy)
    mkpath(pcfg.outdir_xi)

    results, config = load_modal_results(pcfg.input_file)
    results = filter_results(results, pcfg.filter_params, pcfg.filter_deltas)
    grouped = group_results_by_param(results)

    println("Processing entropy and localization parameter ξ with moving average...\n")

    for param_val in sort(collect(keys(grouped)))
        plot_entropy_xi_for_param(
            param_val, grouped[param_val], config,
            pcfg.smooth_delta, pcfg.outdir_entropy, pcfg.outdir_xi)
    end

    println("✓ Done. Entropy and ξ plots saved to:")
    println("  • $(pcfg.outdir_entropy)")
    println("  • $(pcfg.outdir_xi)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
