using Plots, LaTeXStrings

# Definición de funciones
omega_plus(k, delta) = begin
    A = sqrt(1 - (1 - delta^2) * sin(k/2)^2) 
    sqrt(2 + 2A)
end

omega_minus(k, delta) = begin
    A = sqrt(1 - (1 - delta^2) * sin(k/2)^2)
    sqrt(2 - 2A)
end

# Delta para la condición de resonancia
Delta(k1, k2, delta) = begin
    k3 = -k1 - k2
    omega_plus(k3, delta) - omega_minus(k2, delta) - omega_minus(k1, delta)
end

# parámetros
N = 256
delta = 0.4

# crear malla
k1_vals = range(0, 2*π, length=N)
k2_vals = range(0, 2*π, length=N)   

# calcular Delta en la malla
function compute_D(k1_vals, k2_vals, delta)
    Delta.(k1_vals, k2_vals', delta)
end

# funciones de trazado
function plot_D!(p, N, delta; label=nothing, col=:red, show_colorbar=false, ls=:solid, xmax=2*π,
                 guidefs=16, tickfs=12, legendfs=12)
    k1_vals = range(0, xmax, length=N)
    k2_vals = range(0, xmax, length=N)
    D = compute_D(k1_vals, k2_vals, delta)

    # Compute ticks dynamically based on xmax
    all_ticks = [0,  π/2, π, 3π/2, 2π]
    all_labels = [L"0", L"\pi/2", L"\pi", L"3\pi/2", L"2\pi"]
    
    # Filter ticks that are within [0, xmax]
    tick_indices = findall(t -> t <= xmax, all_ticks)
    ticks = (all_ticks[tick_indices], all_labels[tick_indices])

    contour!(p, k1_vals, k2_vals, D;
             levels=[0],
             color=col,
             linewidth=2,
             linestyle=ls,
             xlabel=L"k_{2}",
             ylabel=L"k_{1}",
             xticks=ticks,
             yticks=ticks,
             xlim=(0, xmax),
             ylim=(0, xmax),
             grid=false,
             label=label,
             colorbar=show_colorbar,
             guidefont=font(guidefs),
             tickfont=font(tickfs),
             legendfont=font(legendfs))

    return p
end

function plot_D(N, delta)
    p = plot()
    plot_D!(p, N, delta)
    return p
end

function plot_multiple_D(N, delta_values; colors=[:red, :blue, :green, :orange, :purple],
                         linestyles=[:solid, :dash, :dashdot, :dot, :dashdotdot],
                         xmax=2*π, guidefs=14, tickfs=11, legendfs=12)
    p = plot(legend=:topright, guidefont=font(guidefs), tickfont=font(tickfs), legendfont=font(legendfs))
    for (i, delta) in enumerate(delta_values)
        col = colors[mod1(i, length(colors))]
        ls  = linestyles[mod1(i, length(linestyles))]

        # Dibujamos el contorno sin etiqueta
        plot_D!(p, N, delta; label=nothing, col=col, show_colorbar=(i==1), ls=ls,
                xmax=xmax, guidefs=guidefs, tickfs=tickfs, legendfs=legendfs)

        # Usamos [NaN] con marker=:none para forzar una muestra clara en la leyenda
        plot!(p, [NaN], [NaN], color=col, linewidth=1.4, linestyle=ls, marker=:none,
              label=latexstring("\\Delta \\kappa = $delta"))
    end
    return p
end

# ejecución principal cuando el fichero se corre como script
if abspath(PROGRAM_FILE) == @__FILE__
    plot_multiple_D(N, [0.0, 0.01, 0.05, 0.1])
        mkpath("results/figures/resonance_level_curves")

    savefig("results/figures/resonance_level_curves/level_sets_resonance.pdf")
end
