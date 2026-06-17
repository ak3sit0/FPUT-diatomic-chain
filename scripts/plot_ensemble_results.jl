"""
    plot_ensemble_results.jl

Visualización de resultados del ensamble de fases:
  - Heatmap de energía modal vs tiempo para Δκ ∈ {0.1, 0.3, 0.5}
  - Entropía S̄(t) y dispersión σ(S) vs Δκ
  - Energía óptica relativa E_opt/E_total vs Δκ

Usage:
  julia --project=. scripts/plot_ensemble_results.jl results/data/ensemble_production/ensemble_results_YYYY-MM-DD.jld2
"""

using JLD2, Plots, Statistics, LaTeXStrings
gr()

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

    # --- Figura 1: Heatmaps de energía modal ---
    fig1 = @layout [grid(1,3)]
    p1 = plot(layout=fig1, size=(1500, 400))

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

                hm = heatmap!(p1[i], t, 1:size(modal_E,1), log10.(modal_E .+ 1e-6),
                             title="Δκ = $δ", xlabel="Tiempo (ciclos)", ylabel="Modo",
                             clim=(-6, 1), colorbar=true, cpalette=:viridis)
            end
        end
    end

    savefig(p1, joinpath(outdir, "01_modal_energy_heatmaps.pdf"))
    println("Guardado: $(outdir)/01_modal_energy_heatmaps.pdf")

    # --- Figura 2: Entropía vs tiempo para cada Δκ representativo ---
    fig2 = @layout [grid(1,3)]
    p2 = plot(layout=fig2, size=(1500, 400))

    for (i, δ) in enumerate(delta_targets)
        if haskey(cases_by_delta, δ)
            result = cases_by_delta[δ]
            t = result.scaled_t
            S_mean = result.entropy_mean
            S_std = result.entropy_std

            if length(t) > 0 && length(S_mean) > 0
                plot!(p2[i], t, S_mean, ribbon=S_std, label="S̄(t) ± σ",
                     title="Δκ = $δ", xlabel="Tiempo (ciclos)", ylabel="Entropía",
                     legend=:bottomright, linewidth=2, fillalpha=0.3)
            end
        end
    end

    savefig(p2, joinpath(outdir, "02_entropy_vs_time.pdf"))
    println("Guardado: $(outdir)/02_entropy_vs_time.pdf")

    # --- Figura 3: Resumen global ---
    fig3 = @layout [grid(2,1)]
    p3 = plot(layout=fig3, size=(1000, 800))

    # Panel 3a: Entropía final vs Δκ
    plot!(p3[1], all_deltas, S_final_means, ribbon=S_final_stds,
         title="Entropía final vs Δκ", xlabel="Δκ", ylabel="S̄(∞)",
         legend=:bottomright, linewidth=2, marker=:circle, markersize=6, fillalpha=0.3,
         ylim=(3.5, 4.2))

    # Panel 3b: Energía óptica relativa vs Δκ
    plot!(p3[2], all_deltas, E_opt_finals,
         title="Fracción de energía óptica vs Δκ", xlabel="Δκ", ylabel="E_opt / E_total",
         legend=false, linewidth=2, marker=:circle, markersize=6,
         ylim=(0.2, 0.55), hline=[0.5], line=(:dash, :gray))

    savefig(p3, joinpath(outdir, "03_global_summary.pdf"))
    println("Guardado: $(outdir)/03_global_summary.pdf")

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
