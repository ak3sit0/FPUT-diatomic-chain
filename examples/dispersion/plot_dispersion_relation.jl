using LinearAlgebra, Plots, LaTeXStrings
gr()

# estética común
default(titlefont=font(16, "Times"), guidefont=font(18, "Times"),
        tickfont=font(16, "Times"), legendfont=font(16, "Times"))

# parámetros del sistema
a = 1.0                                    # lattice constant
m = 1.0                                    # atomic mass

# k‑rangos
k_reduced = range(0.0, stop=π/a, length=1200)

# Δκ para la dispersión (spring constant disorder)
# κ₁ = 1 + Δκ, κ₂ = 1 - Δκ
Δκ_values_dispersion = [0.1, 0.3, 0.6, 0.9]
linestyles = [:solid, :dash, :dot, :dashdot]
colors = [:darkblue, :steelblue, :royalblue, :cornflowerblue]
linewidths = [2.4, 2.9, 3.3, 3.6]   # shades of blue with increasing width

# funciones auxiliares para diatomic chain con spring disorder
function compute_frequencies(Δκ::Float64, k_reduced::AbstractVector,
                             a::Float64, m::Float64)
    # Spring constants: κ₁ = 1 + Δκ, κ₂ = 1 - Δκ
    κ₁ = 1.0 + Δκ
    κ₂ = 1.0 - Δκ
    
    # Diatomic chain dispersion with alternating springs
    # ω² = (κ₁ + κ₂)/m ± √[(κ₁ + κ₂)² - 4κ₁κ₂sin²(ka/2)]/m
    κ_sum = κ₁ + κ₂
    κ_prod = κ₁ * κ₂
    
    discriminant = κ_sum^2 .- 4 .* κ_prod .* sin.(k_reduced .* a ./ 2).^2
    sqrtterm = sqrt.(max.(discriminant, 0.0))
    
    ω2_optical = (κ_sum .+ sqrtterm) ./ m      # optical branch
    ω2_acoustic = (κ_sum .- sqrtterm) ./ m     # acoustic branch
    
    sqrt.(ω2_optical), sqrt.(max.(ω2_acoustic, 0.0))
end

function plot_dispersion(Δκ_values::Vector{Float64},
                         k_reduced::AbstractVector,
                         a::Float64, m::Float64)
    freq_data = map(Δκ -> compute_frequencies(Δκ, k_reduced, a, m),
                    Δκ_values)

    p = plot(xlabel=L"k", ylabel=L"\omega", legend=:left, grid=true, gridalpha=0.3,
             size=(800, 600), margin=5Plots.mm)
    xlims!(p, 0, π/a)
    xticks!(p, [0, π/a/2, π/a],
           [L"0", L"\frac{\pi}{2a}", L"\frac{\pi}{a}"])

    max_ω = maximum(vcat([maximum(ω_opt) for (ω_opt, _) in freq_data]...,
                     [maximum(ω_ac)  for (_, ω_ac) in freq_data]...))
    ylims!(p, 0, max_ω * 1.05)   # +5% headroom
    
    # Plot data
    for (i, (ω_opt, ω_ac)) in enumerate(freq_data)
        ls = linestyles[mod1(i, length(linestyles))]
        col = colors[mod1(i, length(colors))]
        lw = linewidths[mod1(i, length(linewidths))]
        κ₁ = 1.0 + Δκ_values[i]
        κ₂ = 1.0 - Δκ_values[i]
        #lbl = L"\Delta \kappa = %$(Δκ_values[i]) \quad (\kappa_1, \kappa_2) = (%.1f, %.1f)" |> 
        #      x -> replace(x, "%.1f" => string(round(κ₁, digits=1)))
        
        # Legend entry with invisible line
        plot!(p, [0, 0], [0, 0], label=L"\Delta \kappa = %$(Δκ_values[i])", 
              linestyle=ls, color=col, linewidth=1.2, marker=:none)

        # Actual curves without label
        plot!(p, k_reduced, ω_opt, label="", linestyle=ls, color=col, 
              linewidth=lw, marker=:none)
        plot!(p, k_reduced, ω_ac, label="", linestyle=ls, color=col, 
              linewidth=lw, marker=:none)
    end

    mkpath("results/figures/dispersion_relation")
    savefig(p, "results/figures/dispersion_relation/dispersion_relation_delta_k.pdf")
    return p
end

# ejecución cuando se llama como script
if abspath(PROGRAM_FILE) == @__FILE__
    plot_dispersion(Δκ_values_dispersion, k_reduced, a, m)
end