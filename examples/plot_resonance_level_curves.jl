using Plots, LaTeXStrings

# Dispersión — forma compacta ec. (9) del paper
# κ* = 1 - Δκ², válido para cualquier convención A/B
function omega_plus(k, delta)
    k = mod(k + π, 2π) - π          # wrap a [-π, π]
    A = sqrt(max(1 - (1 - delta^2) * sin(k/2)^2, 0.0))
    sqrt(2 + 2A)
end

function omega_minus(k, delta)
    k = mod(k + π, 2π) - π
    A = sqrt(max(1 - (1 - delta^2) * sin(k/2)^2, 0.0))
    sqrt(max(2 - 2A, 0.0))
end

# Residuo de la condición de resonancia ω₋(k₁) + ω₋(k₂) = ω₊(k₃)
function resonance_residual(k1, k2, delta)
    k3 = mod(-k1 - k2 + π, 2π) - π  # wrap correcto
    omega_plus(k3, delta) - omega_minus(k1, delta) - omega_minus(k2, delta)
end
 
function compute_residual(k1_vals, k2_vals, delta)
    # Si delta es exactamente 0, forzamos un residuo que NUNCA sea 0
    # para que las curvas de nivel [0.0] salgan completamente vacías.
    if delta == 0.0
        # Devolvemos una matriz llena de un valor constante (ej. 1.0)
        # Así contour! no encontrará ningún cero y la gráfica quedará limpia.
        return fill(1.0, length(k1_vals), length(k2_vals))
    else
        # Para cualquier otro caso (delta > 0), el código sigue igual que antes
        return resonance_residual.(k1_vals, k2_vals', delta)
    end
end


function plot_resonance!(p, N, delta;
                         label=nothing, col=:red, ls=:solid,
                         xmin=-π, xmax=π,
                         guidefs=14, tickfs=11, legendfs=12)
    k1_vals = range(xmin, xmax, length=N)
    k2_vals = range(xmin, xmax, length=N)
    D = compute_residual(k1_vals, k2_vals, delta)

    all_ticks  = [-π, -π/2, 0, π/2, π]
    all_labels = [L"-\pi", L"-\pi/2", L"0", L"\pi/2", L"\pi"]
    idx = findall(t -> xmin <= t <= xmax, all_ticks)

    contour!(p, k1_vals, k2_vals, D;
             levels      = [0.0],
             color       = col,
             linewidth   = 2,
             linestyle   = ls,
             label       = label,
             xlabel      = L"k_1",
             ylabel      = L"k_2",
             xticks      = (all_ticks[idx], all_labels[idx]),
             yticks      = (all_ticks[idx], all_labels[idx]),
             xlim        = (xmin, xmax),
             ylim        = (xmin, xmax),
             grid        = true,
             colorbar    = false,
             guidefont   = font(guidefs),
             tickfont    = font(tickfs),
             legendfont  = font(legendfs))
    return p
end

function plot_multiple_resonance(N, delta_values;
                                  linestyles = [:solid, :dash, :dashdot, :dot, :dashdotdot],
                                  guidefs=14, tickfs=11, legendfs=12)
    colors = [:darkblue, :steelblue, :cornflowerblue, :deepskyblue, :lightblue]
    p = plot(legend=:topright,
             guidefont  = font(guidefs),
             tickfont   = font(tickfs),
             legendfont = font(legendfs))

    for (i, delta) in enumerate(delta_values)
        col = colors[mod1(i, length(colors))]
        ls  = linestyles[mod1(i, length(linestyles))]

        plot_resonance!(p, N, delta;
                        label    = nothing,
                        col      = col,
                        ls       = ls,
                        guidefs  = guidefs,
                        tickfs   = tickfs,
                        legendfs = legendfs)

        # entrada de leyenda limpia
        plot!(p, [NaN], [NaN],
              color     = col,
              linewidth = 1.4,
              linestyle = ls,
              marker    = :none,
              label     = latexstring("\\Delta\\kappa = $delta"))
    end
    return p
end

if abspath(PROGRAM_FILE) == @__FILE__
    p = plot_multiple_resonance(256, [0.1, 0.3, 0.4, 0.49])
    mkpath("results/figures/resonance_level_curves")
    savefig(p, "results/figures/resonance_level_curves/level_sets_resonance.png")
end