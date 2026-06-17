"""
    plot_ensemble_results.jl

Visualización de resultados del ensamble de fases (estilo CairoMakie):
  - Heatmap de energía modal vs tiempo para Δκ ∈ {0.1, 0.3, 0.5}
  - Entropía S̄(t) y dispersión σ(S) vs tiempo para cada Δκ
  - Resumen global: S̄(∞) y E_opt/E_total vs Δκ

Usage:
  julia --project=. scripts/plot_ensemble_results.jl results/data/ensemble_production/ensemble_results_YYYY-MM-DD.jld2
"""

using JLD2, CairoMakie, Statistics, LaTeXStrings, Colors

# Gradiente custom de blues (rescatado de plot_heatmap_grid.jl)
const BLUES_GRADIENT = cgrad([:white, "#B2D9FF", "#5999F2", "#3359CC", "#0D4CB3"])

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

    for result in results
        push!(all_deltas, result.Delta)
        push!(S_final_means, result.entropy_mean[end])
        push!(S_final_stds, result.entropy_std[end])

        E_total = sum(result.E_acoustic_mean[end]) + sum(result.E_optical_mean[end])
        E_opt_frac = sum(result.E_optical_mean[end]) / E_total
        push!(E_opt_finals, E_opt_frac)
    end

    # Ordenar por Δκ
    idx_sort = sortperm(all_deltas)
    all_deltas = all_deltas[idx_sort]
    S_final_means = S_final_means[idx_sort]
    S_final_stds = S_final_stds[idx_sort]
    E_opt_finals = E_opt_finals[idx_sort]

    # --- Figura 1: Heatmaps de energía modal (CairoMakie) ---
    fig1 = Figure(size=(1500, 450))
    ga1 = fig1[1, 1] = GridLayout()

    for (i, δ) in enumerate(delta_targets)
        if haskey(cases_by_delta, δ)
            result = cases_by_delta[δ]
            modal_E = result.modal_E_mean
            t = result.scaled_t

            if !isempty(modal_E) && length(t) > 0
                # Downsample si es muy grande
                nt = size(modal_E, 2)
                if nt > 1000
                    idx_ds = round.(Int, range(1, nt, length=1000))
                    modal_E = modal_E[:, idx_ds]
                    t = t[idx_ds]
                end

                # Normalizar energía a [0, 1] por cada Δκ
                E_min = minimum(modal_E)
                E_max = maximum(modal_E)
                E_norm = (modal_E .- E_min) ./ (E_max - E_min + 1e-12)

                ax = Axis(ga1[1, i];
                    title=latexstring("\\Delta\\kappa = $δ"),
                    xlabel="Tiempo (ciclos)", ylabel="Modo",
                    titlesize=16, xlabelsize=12, ylabelsize=12)

                hm = heatmap!(ax, t, 1:size(E_norm,1), E_norm';
                    colormap=BLUES_GRADIENT, colorrange=(0, 1))

                if i == length(delta_targets)
                    Colorbar(fig1[1, i+1], hm; label="Energía normalizada")
                end
            end
        end
    end

    save(joinpath(outdir, "01_modal_energy_heatmaps.png"), fig1; px_per_unit=2)
    println("Guardado: $(outdir)/01_modal_energy_heatmaps.png")

    # --- Figura 2: Entropía vs tiempo para cada Δκ representativo ---
    fig2 = Figure(size=(1500, 450))
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
                    xlabel="Tiempo (ciclos)", ylabel="Entropía espectral",
                    titlesize=16, xlabelsize=12, ylabelsize=12)

                band!(ax, t, S_mean .- S_std, S_mean .+ S_std, alpha=0.3, color="#5999F2")
                lines!(ax, t, S_mean, color="#0D4CB3", linewidth=2)
            end
        end
    end

    save(joinpath(outdir, "02_entropy_vs_time.png"), fig2; px_per_unit=2)
    println("Guardado: $(outdir)/02_entropy_vs_time.png")

    # --- Figura 3: Resumen global (2 paneles) ---
    fig3 = Figure(size=(1000, 800))
    ga3 = fig3[1, 1] = GridLayout()

    # Panel 3a: Entropía final vs Δκ
    ax_s = Axis(ga3[1, 1];
        title="Entropía final vs Δκ",
        xlabel=latexstring("\\Delta\\kappa"), ylabel=latexstring("\\bar{S}(\\infty)"),
        titlesize=14, xlabelsize=12, ylabelsize=12)

    band!(ax_s, all_deltas, S_final_means .- S_final_stds, S_final_means .+ S_final_stds,
        alpha=0.3, color="#5999F2")
    scatter!(ax_s, all_deltas, S_final_means, color="#0D4CB3", markersize=8)
    lines!(ax_s, all_deltas, S_final_means, color="#0D4CB3", linewidth=2)
    ylims!(ax_s, 3.5, 4.2)

    # Panel 3b: Energía óptica relativa vs Δκ
    ax_e = Axis(ga3[2, 1];
        title="Fracción de energía óptica vs Δκ",
        xlabel=latexstring("\\Delta\\kappa"), ylabel=latexstring("E_{\\mathrm{opt}} / E_{\\mathrm{tot}}"),
        titlesize=14, xlabelsize=12, ylabelsize=12)

    scatter!(ax_e, all_deltas, E_opt_finals, color="#0D4CB3", markersize=8)
    lines!(ax_e, all_deltas, E_opt_finals, color="#0D4CB3", linewidth=2)
    hlines!(ax_e, [0.5], color=:gray, linestyle=:dash, linewidth=1)
    ylims!(ax_e, 0.2, 0.55)

    rowgap!(ga3, 15)
    save(joinpath(outdir, "03_global_summary.png"), fig3; px_per_unit=2)
    println("Guardado: $(outdir)/03_global_summary.png")

    println("\n=== Resumen de resultados ===")
    println("Δκ\tS̄(∞)\t\tE_opt/E_tot")
    println("---\t-------\t\t-----------")
    for (δ, S, E_opt) in zip(all_deltas, S_final_means, E_opt_finals)
        println("$(round(δ; digits=2))\t$(round(S; digits=3))\t\t$(round(E_opt; digits=3))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: julia scripts/plot_ensemble_results.jl <jld2_file>")
    plot_ensemble_diagnostics(ARGS[1])
end
