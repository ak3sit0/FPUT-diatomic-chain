"""
    plot_thermalization_time.jl

Visualiza T_th vs Δκ con efecto halo (varianza como gradiente de color).
Calcula σ(T_th) desde E_optical de cada realización si están disponibles.
Estilo: Plots.jl + scatter elegante, fuentes grandes.

Usage:
  julia --project=. examples/plot_thermalization_time.jl results/data/ensemble_production_100real/ensemble_results_YYYY-MM-DD.jld2
"""

using JLD2, Plots, LaTeXStrings, Statistics

gr()

function apply_recovery_style!()
    default(titlefont=font(16), guidefont=font(14), tickfont=font(11), legendfont=font(12))
end

function compute_thermalization_time(t::Vector, E_opt::Vector; threshold=0.9)
    """Tiempo en que E_opt alcanza el 90% de su valor asintótico."""
    if isempty(E_opt) || length(t) != length(E_opt)
        return NaN
    end
    E_inf = mean(E_opt[max(1, round(Int, 0.9*length(E_opt))):end])
    idx = findfirst(e -> e >= threshold * E_inf, E_opt)
    return isnothing(idx) ? NaN : t[idx]
end

function plot_halo_scatter!(p, deltas, means, stds, color_main, color_halo; n_layers=5)
    """
    Dibuja scatter con efecto de halo (varianza como capas de color degradadas).
    El tamaño del halo es proporcional a la desviación estándar.
    Parámetros:
      - means, stds: centroides y desviaciones estándar
      - n_layers: número de capas concéntricas del halo
    """
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
                t_th = compute_thermalization_time(result.scaled_t, E_opt_real)
                push!(T_th_realizations, t_th)
            end
            valid_th = .!isnan.(T_th_realizations)
            if any(valid_th)
                T_th_mean = mean(T_th_realizations[valid_th])
                T_th_std = std(T_th_realizations[valid_th])
            end
        else
            # Fallback: solo media del ensamble (compatibilidad con datos antiguos)
            T_th_mean = compute_thermalization_time(result.scaled_t, result.E_optical_mean)
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

    # --- Figura 1: T_th vs Δκ (escala log) con halo ---
    p1 = plot(
        xlabel=latexstring("\\Delta\\kappa"),
        ylabel=latexstring("T_{\\mathrm{th}} \\; (\\mathrm{cycles})"),
        legend=false,
        grid=true,
        size=(1000, 700),
        yscale=:log10,
        framestyle=:box,
        bottom_margin=5Plots.mm,
        left_margin=5Plots.mm
    )

    plot_halo_scatter!(p1, all_deltas, T_th_means, T_th_stds,
        :darkblue, :steelblue; n_layers=5)

    savefig(p1, "results/figures/ensemble/thermalization_time_log.pdf")
    println("Saved: results/figures/ensemble/thermalization_time_log.pdf")

    # --- Figura 2: T_th vs Δκ (escala lineal) con halo ---
    p2 = plot(
        xlabel=latexstring("\\Delta\\kappa"),
        ylabel=latexstring("T_{\\mathrm{th}} \\; (\\mathrm{cycles})"),
        legend=false,
        grid=true,
        size=(1000, 700),
        framestyle=:box,
        bottom_margin=5Plots.mm,
        left_margin=5Plots.mm
    )

    plot_halo_scatter!(p2, all_deltas, T_th_means, T_th_stds,
        :steelblue, :lightblue; n_layers=5)

    savefig(p2, "results/figures/ensemble/thermalization_time_linear.pdf")
    println("Saved: results/figures/ensemble/thermalization_time_linear.pdf")

    # --- Figura 3: E_opt equilibrium vs Δκ con halo ---
    p3 = plot(
        xlabel=latexstring("\\Delta\\kappa"),
        ylabel=latexstring("E_{\\mathrm{opt}} / E_{\\mathrm{tot}}"),
        legend=false,
        grid=true,
        size=(1000, 700),
        framestyle=:box,
        ylim=(0.2, 0.6),
        bottom_margin=5Plots.mm,
        left_margin=5Plots.mm
    )

    plot_halo_scatter!(p3, all_deltas, E_opt_means, E_opt_stds,
        :darkred, :salmon; n_layers=3)
    hline!(p3, [0.5], line=(:dash, :gray, 2), label="")

    savefig(p3, "results/figures/ensemble/optical_energy_equilibrium.pdf")
    println("Saved: results/figures/ensemble/optical_energy_equilibrium.pdf")

    # --- Figura 4: Ambas en subplots ---
    p_layout = @layout [a; b]

    p4a = plot(
        xlabel="", ylabel=latexstring("T_{\\mathrm{th}}"),
        legend=false, grid=true, yscale=:log10, framestyle=:box,
        bottom_margin=2Plots.mm, left_margin=5Plots.mm
    )
    plot_halo_scatter!(p4a, all_deltas, T_th_means, T_th_stds,
        :darkblue, :steelblue; n_layers=5)

    p4b = plot(
        xlabel=latexstring("\\Delta\\kappa"),
        ylabel=latexstring("E_{\\mathrm{opt}}/E_{\\mathrm{tot}}"),
        legend=false, grid=true, ylim=(0.2, 0.6), framestyle=:box,
        bottom_margin=5Plots.mm, left_margin=5Plots.mm
    )
    plot_halo_scatter!(p4b, all_deltas, E_opt_means, E_opt_stds,
        :darkred, :salmon; n_layers=3)
    hline!(p4b, [0.5], line=(:dash, :gray, 2), label="")

    p_combined = plot(p4a, p4b, layout=p_layout, size=(1000, 1100))
    savefig(p_combined, "results/figures/ensemble/thermalization_combined.pdf")
    println("Saved: results/figures/ensemble/thermalization_combined.pdf")

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
