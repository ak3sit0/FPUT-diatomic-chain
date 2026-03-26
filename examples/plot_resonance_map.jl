using LinearAlgebra, Plots, LaTeXStrings 

default(titlefont=font(14, "Times"), guidefont=font(14, "Times"),
        tickfont=font(14, "Times"), legendfont=font(14, "Times"))

function plot_resonance_map(m::Float64)
    Delta_kappa = range(0.0, 0.99, length=500)

    omega_acu_max = sqrt.(2.0 .- 2.0 .* Delta_kappa)
    omega_opt_min = sqrt.(2.0 .+ 2.0 .* Delta_kappa)
    omega_opt_max = fill(2.0, length(Delta_kappa))

    p = plot(xlabel=L"\Delta \kappa", ylabel=L"\omega",
             legend=:topright, grid=false,
             guidefont=font(14, "Times"), tickfont=font(14, "Times"),
             ylim=(0, 4.0), xlim=(0, 0.99), gridalpha=0.3) 

    plot!(p, Delta_kappa, omega_opt_min, color=:red, linewidth=2,
    #      label=L"$\omega_{min}^{+}$") # Add the plot for the optical band minimum frequency
            label="")  
    
            plot!(p, Delta_kappa, omega_opt_min, fillrange=omega_opt_max,
          color=:red, alpha=0.15, label="Optical Band") # Add the shaded area for the optical band between its minimum and maximum frequencies

    plot!(p, Delta_kappa, omega_acu_max, color=:blue, linewidth=2,
    #      label=L"$\omega_{max}^{-}$") # Add the plot for the acoustic band maximum frequency
            label="")

    plot!(p, Delta_kappa, zeros(length(Delta_kappa)),
          fillrange=omega_acu_max, color=:blue, alpha=0.15, # Add the shaded area for the acoustic band between zero and its maximum frequency
          label="Acoustic Band")

    for n in (2, 3)
        harmonic = n .* omega_acu_max # Calculate the nth harmonic of the acoustic band maximum frequency

        plot!(p, Delta_kappa, harmonic, linestyle=:dash, color=:black,
              alpha=0.5, label="")

        diff = harmonic .- omega_opt_min # Calculate the difference between the nth harmonic and the optical band minimum frequency to find resonance crossings
        indices = findall(x -> x >= 0, diff) # Find indices where the difference is non-negative, indicating a crossing of the nth harmonic with the optical band minimum frequency
        if !isempty(indices)
            crossing_idx = indices[end] # Get the last index where the crossing occurs (the highest Delta_kappa value where the nth harmonic is still above the optical band minimum frequency)
            x_pt = Delta_kappa[crossing_idx]
            y_pt = harmonic[crossing_idx]

            scatter!(p, [x_pt], [y_pt], color=:black, markersize=5,
                     label="") # Mark the crossing point on the plot
            annotate!(p, x_pt + 0.01, y_pt + 0.15, text("$n:1", 10, :black))  # Add annotation for the resonance crossing point, indicating the harmonic ratio (n:1) at that point on the plot
        end
    end

    mkpath("results/figures/resonance_map")
    savefig(p, "results/figures/resonance_map/resonance_map_delta_k.pdf")
    return p
end

if abspath(PROGRAM_FILE) == @__FILE__
    m = 1.0
    plot_resonance_map(m)
end