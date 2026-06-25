"""
    plot_entropy_Nsweep.jl

Grafica la evolución de entropía para múltiples N en una sola gráfica.
Basado en plot_entropy_N_timeseries.jl con estilo profesional.

Usage:
  julia --project=. examples/plot_entropy_Nsweep.jl results/data/nsweep_test/nsweep_results_YYYY-MM-DD.jld2
"""

using JLD2, Plots, LaTeXStrings, Statistics, Printf

gr()

const PLOT_MAX_POINTS = 2000
const LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]

function apply_recovery_style!()
    default(titlefont=font(16), guidefont=font(14), tickfont=font(11), legendfont=font(12))
end

function main()
    if isempty(ARGS)
        println("Usage: julia examples/plot_entropy_Nsweep.jl <results.jld2>")
        return
    end

    data = jldopen(ARGS[1], "r")
    results = data["results"]
    close(data)

    # Ordenar por N
    sort!(results, by=r -> r.N)

    apply_recovery_style!()

    # Paleta de colores
    palette = [:blue, :red, :green, :orange, :purple, :brown, :magenta]

    # Crear gráfica
    p = plot(
        xlabel=L"t \quad (\mathrm{cycles})",
        ylabel=L"S(t)",
        lw=2.5,
        legend=:bottomright,
        grid=false,
        size=(1000, 600),
        xscale=:log10,
        framestyle=:box
    )

    # Plotear cada N
    for (i, result) in enumerate(results)
        N = result.N
        t = result.scaled_t
        S = result.entropy

        # Convertir a Float64 por si acaso
        t = Float64.(t)
        S = Float64.(S)

        # Downsample si hay muchos puntos
        npts = length(t)
        if npts > PLOT_MAX_POINTS
            step = max(1, Int(floor(npts / PLOT_MAX_POINTS)))
            idx = 1:step:npts
            t = t[idx]
            S = S[idx]
        end

        color = palette[mod1(i, length(palette))]
        ls = LINESTYLES[mod1(i, length(LINESTYLES))]

        plot!(p, t, S;
            color=color,
            linestyle=ls,
            linewidth=2.2,
            label=latexstring("N = $(N)"))
    end

    # Guardar PDF
    outdir = "results/figures/nsweep"
    mkpath(outdir)
    savefig(p, joinpath(outdir, "entropy_Nsweep.pdf"))
    println("Saved: $(outdir)/entropy_Nsweep.pdf")

    # PNG para preview
    savefig(p, joinpath(outdir, "entropy_Nsweep.png"))
    println("Saved: $(outdir)/entropy_Nsweep.png")

    # Tabla de resumen
    println("\n=== Entropy Summary ===")
    println("N\tS_init\tS_final\tΔS")
    println("-\t------\t-------\t---")
    for result in results
        N = result.N
        S_init = result.entropy[1]
        S_final = result.entropy[end]
        ΔS = S_final - S_init
        @printf "%d\t%.3f\t%.3f\t%.3f\n" N S_init S_final ΔS
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
