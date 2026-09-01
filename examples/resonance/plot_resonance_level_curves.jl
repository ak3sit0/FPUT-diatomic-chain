using Plots, LaTeXStrings, Colors
include("../../src/dispersion.jl");      using .Dispersion
include("../../src/fput_analysis.jl");   using .FPUTAnalysis
include("../../src/plotting_utils.jl");  using .PlottingUtils

# Dispersion in compact form (κ* = 1 - Δκ²) — eq. (9) of the paper, in Dispersion
const omega_plus  = Dispersion.omega_compact_plus
const omega_minus = Dispersion.omega_compact_minus

# Residual of the resonance condition ω₋(k₁) + ω₋(k₂) = ω₊(k₃)
function resonance_residual(k1, k2, delta)
    k3 = mod(-k1 - k2 + π, 2π) - π  # correct wrap
    omega_plus(k3, delta) - omega_minus(k1, delta) - omega_minus(k2, delta)
end

function compute_residual(k1_vals, k2_vals, delta)
    # If delta is exactly 0, force a residual that is NEVER 0 so the [0.0]
    # level curves come out completely empty.
    if delta == 0.0
        # Return a matrix filled with a constant value (e.g. 1.0), so
        # contour! finds no zero and the plot stays clean.
        return fill(1.0, length(k1_vals), length(k2_vals))
    else
        # For any other case (delta > 0), the code proceeds as before.
        return resonance_residual.(k1_vals, k2_vals', delta)
    end
end

function add_klapp_umklapp_regions!(p; xmin=-π, xmax=π)
    tri_topright = Shape([0.0, π, π], [π, π, 0.0])
    tri_botleft  = Shape([0.0, -π, -π], [-π, -π, 0.0])

    plot!(p, tri_topright; fillcolor=RGB(0.961, 0.773, 0.094), fillalpha=0.18,
          linealpha=0, label="")
    plot!(p, tri_botleft;  fillcolor=RGB(0.961, 0.773, 0.094), fillalpha=0.18,
          linealpha=0, label="")

    # Boundary k1+k2 = ±π
    plot!(p, [0.0, π], [π, 0.0]; color=RGB(0.753, 0.439, 0.0), linestyle=:dash,
          linewidth=2.5, alpha=0.90, label="")
    plot!(p, [0.0, -π], [-π, 0.0]; color=RGB(0.753, 0.439, 0.0), linestyle=:dash,
          linewidth=2.5, alpha=0.90, label="")

    annotate!(p, (-2.95, -2.35, text("Umklapp\n" * L"|k_1+k_2|>\pi", 9,
              :firebrick, :left)))
    annotate!(p, (-1.15, 2.35, text("Klapp (normal)\n" * L"|k_1+k_2|<\pi", 9,
              :navy, :left)))
    return p
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
             framestyle = :origin,
             colorbar    = false,
             guidefont   = font(guidefs),
             tickfont    = font(tickfs),
             legendfont  = font(legendfs))
    return p
end

function plot_multiple_resonance(N, delta_values;
                                  linestyles = [:solid, :dash, :dashdot, :dot, :dashdotdot],
                                  guidefs=14, tickfs=11, legendfs=12)
    p = plot(legend=:topright,
             guidefont  = font(guidefs),
             tickfont   = font(tickfs),
             legendfont = font(legendfs),
             framestyle = :origin)

    add_klapp_umklapp_regions!(p)   # ← NEW line, before the loop

    for (i, delta) in enumerate(delta_values)
        col = cyc(PALETTE_DELTA, i)
        ls  = cyc(linestyles, i)

        plot_resonance!(p, N, delta;
                        label    = nothing,
                        col      = col,
                        ls       = ls,
                        guidefs  = guidefs,
                        tickfs   = tickfs,
                        legendfs = legendfs)

        # clean legend entry
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
    p = plot_multiple_resonance(256, [0.05, 0.1, 0.3, 0.49])
    mkpath("results/figures/resonance_level_curves")
    savefig(p, "results/figures/resonance_level_curves/level_sets_resonance.pdf")
end