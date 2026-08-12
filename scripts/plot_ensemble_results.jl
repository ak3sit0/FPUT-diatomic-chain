"""
    plot_ensemble_results.jl

Visualización de resultados del ensamble de fases (estilo CairoMakie, layouts corregidos):
  - Heatmap de energía modal vs tiempo para Δκ ∈ {0.1, 0.3, 0.5}
  - Entropía S̄(t) y dispersión σ(S) vs tiempo para cada Δκ
  - Resumen global: S̄(∞) y E_opt/E_total vs Δκ
  - Tiempo de termalización T_th vs Δκ

Usage:
  julia --project=. scripts/plot_ensemble_results.jl results/data/ensemble_production/ensemble_results_YYYY-MM-DD.jld2
"""

using JLD2, CairoMakie, Statistics, LaTeXStrings, Colors

const BLUES_GRADIENT = cgrad([:white, "#B2D9FF", "#5999F2", "#3359CC", "#0D4CB3"])

function compute_thermalization_time(t::Vector, E_opt::Vector; threshold=0.9)
    """Calcula el tiempo en que E_opt alcanza el 90% de su valor asintótico."""
    if isempty(E_opt) || length(t) != length(E_opt)
        return NaN
    end
    E_inf = mean(E_opt[max(1, round(Int, 0.9*length(E_opt))):end])
    idx = findfirst(e -> e >= threshold * E_inf, E_opt)
    return isnothing(idx) ? NaN : t[idx]
end

function plot_ensemble_diagnostics(jld2_path::String; outdir="results/figures/ensemble")
    mkpath(outdir)

    # Cargar resultados
    data = jldopen(jld2_path, "r")
    results = data["results"]
    close(data)

    # Seleccionar casos representativos
    delta_targets = [0.1, 0.3, 0.5]
    cases_by_delta = Dict()

    for result in results
        δ = result.Delta
        if δ in delta_targets
            if !haskey(cases_by_delta, δ)
                cases_by_delta[δ] = result
            end
        end
    end

    # Recolectar todos los Δκ y sus estadísticas
    all_deltas = Float64[]
    S_final_means = Float64[]
    S_final_stds = Float64[]
    E_opt_finals = Float64[]
    T_th_values = Float64[]

    for result in results
        push!(all_deltas, result.Delta)
        push!(S_final_means, result.entropy_mean[end])
        push!(S_final_stds, result.entropy_std[end])

        E_total = sum(result.E_acoustic_mean[end]) + sum(result.E_optical_mean[end])
        E_opt_frac = sum(result.E_optical_mean[end]) / E_total
        push!(E_opt_finals, E_opt_frac)

        # Calcular tiempo de termalización
        t_th = compute_thermalization_time(result.scaled_t, result.E_optical_mean)
        push!(T_th_values, t_th)
    end

    # Ordenar por Δκ
    idx_sort = sortperm(all_deltas)
    all_deltas = all_deltas[idx_sort]
    S_final_means = S_final_means[idx_sort]
    S_final_stds = S_final_stds[idx_sort]
    E_opt_finals = E_opt_finals[idx_sort]
    T_th_values = T_th_values[idx_sort]

    # --- Figura 1: Heatmaps de energía modal ---
    n_cols = length(delta_targets)
    fig1 = Figure(size=(500*n_cols + 150, 500))
    ga1 = fig1[1, 1] = GridLayout()

    hm_ref = nothing

    for (i, δ) in enumerate(delta_targets)
        if haskey(cases_by_delta, δ)
            result = cases_by_delta[δ]
            modal_E = result.modal_E_mean
            t = result.scaled_t

            if !isempty(modal_E) && length(t) > 0
                nt = size(modal_E, 2)
                if nt > 1000
                    idx_ds = round.(Int, range(1, nt, length=1000))
                    modal_E = modal_E[:, idx_ds]
                    t = t[idx_ds]
                end

                E_min = minimum(modal_E)
                E_max = maximum(modal_E)
                E_norm = (modal_E .- E_min) ./ (E_max - E_min + 1e-12)

                ax = Axis(ga1[1, i];
                    title=latexstring("\\Delta\\kappa = $δ"),
                    xlabel="", ylabel="Mode",
                    titlesize=16, xlabelsize=11, ylabelsize=11,
                    xticklabelsvisible = (i == n_cols),
                    yticklabelsvisible = (i == 1))

                hm = heatmap!(ax, t, 1:size(E_norm,1), E_norm';
                    colormap=BLUES_GRADIENT, colorrange=(0, 1))

                if isnothing(hm_ref)
                    hm_ref = hm
                end
            end
        end
    end

    # Colorbar en columna 2 de fig1 (hermana de ga1)
    Colorbar(fig1[1, 2], hm_ref; label="Normalized energy", width=20, ticklabelsize=11)

    # Etiqueta global eje X
    Label(fig1[2, 1], "Time (cycles)"; fontsize=13)

    # Espaciado entre paneles
    colgap!(ga1, 15)

    save(joinpath(outdir, "01_modal_energy_heatmaps.png"), fig1; px_per_unit=2)
    println("Saved: $(outdir)/01_modal_energy_heatmaps.png")

    # --- Figura 2: Entropía vs tiempo ---
    fig2 = Figure(size=(500*n_cols + 150, 500))
    ga2 = fig2[1, 1] = GridLayout()

    for (i, δ) in enumerate(delta_targets)
        if haskey(cases_by_delta, δ)
            result = cases_by_delta[δ]
            t = result.scaled_t
            S_mean = result.entropy_mean
            S_std = result.entropy_std

            if length(t) > 0 && length(S_mean) > 0
                ax = Axis(ga2[1, i];
                    title=latexstring("\\Delta\\kappa = $δ"),
                    xlabel="", ylabel="Spectral entropy",
                    titlesize=16, xlabelsize=11, ylabelsize=11,
                    xticklabelsvisible = (i == n_cols),
                    yticklabelsvisible = (i == 1))

                band!(ax, t, S_mean .- S_std, S_mean .+ S_std, alpha=0.3, color="#5999F2")
                lines!(ax, t, S_mean, color="#0D4CB3", linewidth=2)
            end
        end
    end

    Label(fig2[2, 1], "Time (cycles)"; fontsize=13)
    colgap!(ga2, 15)

    save(joinpath(outdir, "02_entropy_vs_time.png"), fig2; px_per_unit=2)
    println("Saved: $(outdir)/02_entropy_vs_time.png")

    # --- Figura 3: Resumen global ---
    fig3 = Figure(size=(1000, 800))
    ga3 = fig3[1, 1] = GridLayout()

    # Panel 3a: Entropía final vs Δκ
    ax_s = Axis(ga3[1, 1];
        title="Final entropy vs Δκ",
        xlabel=latexstring("\\Delta\\kappa"), ylabel=latexstring("\\bar{S}(\\infty)"),
        titlesize=14, xlabelsize=12, ylabelsize=12)

    band!(ax_s, all_deltas, S_final_means .- S_final_stds, S_final_means .+ S_final_stds,
        alpha=0.3, color="#5999F2")
    scatter!(ax_s, all_deltas, S_final_means, color="#0D4CB3", markersize=8)
    lines!(ax_s, all_deltas, S_final_means, color="#0D4CB3", linewidth=2)
    ylims!(ax_s, 3.5, 4.2)

    # Panel 3b: Energía óptica relativa vs Δκ
    ax_e = Axis(ga3[2, 1];
        title="Optical energy fraction vs Δκ",
        xlabel=latexstring("\\Delta\\kappa"), ylabel=latexstring("E_{\\mathrm{opt}} / E_{\\mathrm{tot}}"),
        titlesize=14, xlabelsize=12, ylabelsize=12)

    scatter!(ax_e, all_deltas, E_opt_finals, color="#0D4CB3", markersize=8)
    lines!(ax_e, all_deltas, E_opt_finals, color="#0D4CB3", linewidth=2)
    hlines!(ax_e, [0.5], color=:gray, linestyle=:dash, linewidth=1)
    ylims!(ax_e, 0.2, 0.55)

    rowgap!(ga3, 20)
    save(joinpath(outdir, "03_global_summary.png"), fig3; px_per_unit=2)
    println("Saved: $(outdir)/03_global_summary.png")

    # --- Figura 4: Tiempo de termalización ---
    fig4 = Figure(size=(1000, 800))
    ga4 = fig4[1, 1] = GridLayout()

    # Panel 4a: T_th vs Δκ
    ax_th = Axis(ga4[1, 1];
        title="Thermalization time vs Δκ",
        xlabel=latexstring("\\Delta\\kappa"), ylabel=latexstring("T_{\\mathrm{th}} \\; (\\mathrm{cycles})"),
        titlesize=14, xlabelsize=12, ylabelsize=12,
        yscale=log10)

    valid_th = .!isnan.(T_th_values)
    if any(valid_th)
        scatter!(ax_th, all_deltas[valid_th], T_th_values[valid_th], color="#0D4CB3", markersize=8)
        lines!(ax_th, all_deltas[valid_th], T_th_values[valid_th], color="#0D4CB3", linewidth=2)
    end

    # Panel 4b: E_opt equilibrium vs Δκ
    ax_eq = Axis(ga4[2, 1];
        title="Equilibrium optical energy fraction",
        xlabel=latexstring("\\Delta\\kappa"), ylabel=latexstring("E_{\\mathrm{opt}}^{\\mathrm{eq}} / E_{\\mathrm{tot}}"),
        titlesize=14, xlabelsize=12, ylabelsize=12)

    scatter!(ax_eq, all_deltas, E_opt_finals, color="#0D4CB3", markersize=8)
    lines!(ax_eq, all_deltas, E_opt_finals, color="#0D4CB3", linewidth=2)
    hlines!(ax_eq, [0.5], color=:gray, linestyle=:dash, linewidth=1)
    ylims!(ax_eq, 0.2, 0.55)

    rowgap!(ga4, 20)
    save(joinpath(outdir, "04_thermalization_time.png"), fig4; px_per_unit=2)
    println("Saved: $(outdir)/04_thermalization_time.png")

    println("\n=== Thermalization times ===")
    println("Δκ\tT_th (cycles)\tE_opt/E_tot")
    println("---\t-----------\t-----------")
    for (δ, t_th, e_opt) in zip(all_deltas, T_th_values, E_opt_finals)
        t_str = isnan(t_th) ? "∞" : "$(round(t_th; digits=1))"
        println("$(round(δ; digits=2))\t$t_str\t\t$(round(e_opt; digits=3))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: julia scripts/plot_ensemble_results.jl <jld2_file>")
    plot_ensemble_diagnostics(ARGS[1])
end
