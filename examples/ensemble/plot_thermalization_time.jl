"""
    plot_thermalization_time.jl

Visualiza T_th vs Δκ con efecto halo (varianza como gradiente de color).
Calcula σ(T_th) desde E_optical de cada realización si están disponibles.
Estilo: Plots.jl + scatter elegante, fuentes grandes.

Usage:
  julia --project=. examples/plot_thermalization_time.jl results/data/ensemble_production_100real/ensemble_results_YYYY-MM-DD.jld2
"""

using JLD2, Plots, LaTeXStrings, Statistics
include("../../src/fput_analysis.jl");   using .FPUTAnalysis
include("../../src/plotting_utils.jl");  using .PlottingUtils

gr()

const XLABEL_DELTA = latexstring("\\Delta\\kappa")
const YLABEL_TTH   = latexstring("T_{\\mathrm{th}} \\; (\\mathrm{cycles})")
const YLABEL_EOPT  = latexstring("E_{\\mathrm{opt}} / E_{\\mathrm{tot}}")

function apply_recovery_style!()
    default(titlefont=font(16), guidefont=font(14), tickfont=font(11), legendfont=font(12))
end

"""
    plot_halo_scatter!(p, deltas, means, stds, color_main, color_halo; n_layers=5)

Dibuja scatter con efecto de halo (varianza como capas de color degradadas).
El tamaño del halo es proporcional a la desviación estándar.
"""
function plot_halo_scatter!(p, deltas, means, stds, color_main, color_halo; n_layers=5)
    max_std = maximum(stds[.!isnan.(stds)])
    max_std = max_std > 0 ? max_std : 1.0

    for layer in 1:n_layers
        alpha_halo = 0.25 * (1.0 - (layer - 1) / n_layers)
        # Tamaño base 6, más capas proporcionales a σ normalizado
        markersize_layer = 6 .+ (stds ./ max_std) .* layer .* 3

        scatter!(p, deltas, means;
            marker=:circle,
            markersize=markersize_layer,
            markerstrokewidth=0,
            markerstrokecolor=color_halo,
            color=color_halo,
            alpha=alpha_halo,
            label="",
            legend=false)
    end

    scatter!(p, deltas, means;
        marker=:circle,
        markersize=6,
        markerstrokewidth=0,
        markerstrokecolor=color_main,
        color=color_main,
        label="",
        legend=false)
end

"""
    halo_panel(deltas, means, stds; colors, ylabel, n_layers, guide, kwargs...) -> Plot

Panel de scatter con halo. Las cuatro figuras de este script son este mismo panel
con distinta escala, paleta y línea guía.
"""
function halo_panel(deltas, means, stds; colors, ylabel, n_layers=5,
                    guide=nothing, xlabel=XLABEL_DELTA, kwargs...)
    p = plot(; xlabel=xlabel, ylabel=ylabel, legend=false, grid=true,
               framestyle=:box, left_margin=5Plots.mm, kwargs...)
    plot_halo_scatter!(p, deltas, means, stds, colors...; n_layers=n_layers)
    isnothing(guide) || hline!(p, [guide], line=(:dash, GRAY_GUIDE, 2), label="")
    p
end

function main()
    if isempty(ARGS)
        println("Usage: julia examples/plot_thermalization_time.jl <results.jld2>")
        return
    end

    data = jldopen(ARGS[1], "r")
    results = data["results"]
    close(data)

    # Recolectar datos por Δκ
    all_deltas = Float64[]
    T_th_means = Float64[]
    T_th_stds = Float64[]
    E_opt_means = Float64[]
    E_opt_stds = Float64[]

    for result in results
        push!(all_deltas, result.Delta)

        # Calcular T_th: intentar con per-realización, fallback a media
        T_th_mean = NaN
        T_th_std = NaN

        if haskey(result, :E_optical_realizations) && !isempty(result.E_optical_realizations)
            # Tenemos datos per-realización: calcular σ(T_th)
            T_th_realizations = Float64[]
            for E_opt_real in result.E_optical_realizations
                t_th = thermalization_time(result.scaled_t, E_opt_real)
                push!(T_th_realizations, t_th)
            end
            valid_th = .!isnan.(T_th_realizations)
            if any(valid_th)
                T_th_mean = mean(T_th_realizations[valid_th])
                T_th_std = std(T_th_realizations[valid_th])
            end
        else
            # Fallback: solo media del ensamble (compatibilidad con datos antiguos)
            T_th_mean = thermalization_time(result.scaled_t, result.E_optical_mean)
            T_th_std = 0.0
        end

        push!(T_th_means, T_th_mean)
        push!(T_th_stds, T_th_std)

        # E_opt equilibrio (media del ensamble)
        E_total = sum(result.E_acoustic_mean[end]) + sum(result.E_optical_mean[end])
        E_opt_frac = sum(result.E_optical_mean[end]) / E_total
        push!(E_opt_means, E_opt_frac)
        push!(E_opt_stds, 0.0)
    end

    # Ordenar por Δκ
    idx_sort = sortperm(all_deltas)
    all_deltas = all_deltas[idx_sort]
    T_th_means = T_th_means[idx_sort]
    T_th_stds = T_th_stds[idx_sort]
    E_opt_means = E_opt_means[idx_sort]
    E_opt_stds = E_opt_stds[idx_sort]

    apply_recovery_style!()

    outdir = "results/figures/ensemble"
    mkpath(outdir)

    T_th_panel(colors; kwargs...) = halo_panel(all_deltas, T_th_means, T_th_stds;
        colors=colors, ylabel=YLABEL_TTH, n_layers=5, kwargs...)
    E_opt_panel(; kwargs...) = halo_panel(all_deltas, E_opt_means, E_opt_stds;
        colors=HALO_RED, ylabel=YLABEL_EOPT, n_layers=3, guide=0.5, ylim=(0.2, 0.6), kwargs...)

    figs = [
        ("thermalization_time_log.pdf",
         T_th_panel(HALO_BLUE; yscale=:log10, size=(1000, 700), bottom_margin=5Plots.mm)),
        ("thermalization_time_linear.pdf",
         T_th_panel(HALO_TEAL; size=(1000, 700), bottom_margin=5Plots.mm)),
        ("optical_energy_equilibrium.pdf",
         E_opt_panel(size=(1000, 700), bottom_margin=5Plots.mm)),
        # Figura 4: las dos anteriores apiladas (sin xlabel en el panel superior)
        ("thermalization_combined.pdf",
         plot(T_th_panel(HALO_BLUE; yscale=:log10, xlabel="",
                         ylabel=latexstring("T_{\\mathrm{th}}"), bottom_margin=2Plots.mm),
              E_opt_panel(bottom_margin=5Plots.mm),
              layout=(@layout [a; b]), size=(1000, 1100))),
    ]

    for (name, fig) in figs
        savefig(fig, joinpath(outdir, name))
        println("Saved: $outdir/$name")
    end

    # --- Tabla de resultados ---
    println("\n=== Thermalization times with uncertainty ===")
    println("Δκ\t\tT_th ± σ(T_th) [cycles]\t\tE_opt/E_tot")
    println("---\t\t------------------------\t\t-----------")
    for (δ, t_mean, t_std, e_opt) in zip(all_deltas, T_th_means, T_th_stds, E_opt_means)
        if isnan(t_mean)
            t_str = "∞"
        else
            t_str = "$(round(t_mean; digits=1)) ± $(round(t_std; digits=1))"
        end
        println("$(round(δ; digits=2))\t\t$t_str\t$(round(e_opt; digits=3))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
