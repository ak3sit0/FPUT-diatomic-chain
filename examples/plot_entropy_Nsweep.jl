"""
    plot_entropy_Nsweep.jl

Grafica la evolución de entropía para múltiples N en una sola gráfica.
Estilo: Plots.jl + gr(), curvas superpuestas, colores y estilos distintos por N.

Usage:
  julia --project=. examples/plot_entropy_Nsweep.jl results/data/nsweep_test/nsweep_results_YYYY-MM-DD.jld2
"""

using JLD2, Plots, LaTeXStrings, Statistics

gr()

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

    # Colores y estilos por N
    palette = [:darkblue, :darkred, :darkgreen, :darkorange]
    linestyles = [:solid, :dash, :dot, :dashdot]
    markers = [:circle, :square, :diamond, :pentagon]

    # Crear gráfica
    p = plot(
        xlabel=latexstring("t \\; (\\mathrm{cycles})"),
        ylabel=latexstring("S(t)"),
        title="Spectral entropy vs time for different N",
        legend=:bottomright,
        lw=2.5,
        grid=true,
        size=(1000, 700),
        xscale=:log10,
        framestyle=:box,
        bottom_margin=5Plots.mm,
        left_margin=5Plots.mm
    )

    # Plotear cada N
    for (i, result) in enumerate(results)
        N = result.N
        t = result.scaled_t
        S = result.entropy

        # Downsample si hay muchos puntos
        npts = length(t)
        if npts > 2000
            step = max(1, Int(floor(npts / 2000)))
            idx = 1:step:npts
            t = t[idx]
            S = S[idx]
        end

        color = palette[mod1(i, length(palette))]
        ls = linestyles[mod1(i, length(linestyles))]
        marker = markers[mod1(i, length(markers))]

        plot!(p, t, S;
            color=color,
            linestyle=ls,
            marker=marker,
            markersize=4,
            markerstrokewidth=0,
            lw=2.5,
            label=latexstring("N = $(N)"),
            legend=:bottomright)
    end

    # Línea horizontal en log(N) para referencia
    N_values = [r.N for r in results]
    if !isempty(N_values)
        S_max = maximum(N_values) |> log
        hline!(p, [S_max], line=(:dash, :gray, 1.5), label=latexstring("\\log N_{\\max}"))
    end

    # Guardar
    outdir = "results/figures/nsweep"
    mkpath(outdir)
    savefig(p, joinpath(outdir, "entropy_Nsweep.pdf"))
    println("Saved: $(outdir)/entropy_Nsweep.pdf")

    # Tabla de resumen
    println("\n=== Entropy Summary ===")
    println("N\tS_init\tS_final\tΔS")
    println("-\t------\t-------\t---")
    for result in results
        N = result.N
        S_init = result.entropy[1]
        S_final = result.entropy[end]
        ΔS = S_final - S_init
        println("$(N)\t$(round(S_init; digits=3))\t$(round(S_final; digits=3))\t$(round(ΔS; digits=3))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
