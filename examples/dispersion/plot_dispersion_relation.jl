using LinearAlgebra, Plots, LaTeXStrings
include("../../src/dispersion.jl");      using .Dispersion
include("../../src/fput_analysis.jl");   using .FPUTAnalysis
include("../../src/plotting_utils.jl");  using .PlottingUtils
gr()

# Common style
default(titlefont=font(16, "Times"), guidefont=font(18, "Times"),
        tickfont=font(16, "Times"), legendfont=font(16, "Times"))

# System parameters
a = 1.0                                    # lattice constant
m = 1.0                                    # atomic mass

# k‑rangos
k_reduced = range(0.0, stop=π/a, length=1200)

# Δκ for the dispersion (spring constant disorder)
# κ₁ = 1 + Δκ, κ₂ = 1 - Δκ
Δκ_values_dispersion = [0.1, 0.3, 0.6, 0.9]

"""
    compute_frequencies(Δκ, k, a, m) -> (ω_optical, ω_acoustic)

Ambas ramas para κ₁ = 1 + Δκ, κ₂ = 1 - Δκ. La dispersión vive en
`src/dispersion.jl`; aquí sólo se reescala por la masa (m = 1 en `Dispersion`).
"""
function compute_frequencies(Δκ::Float64, k_reduced::AbstractVector,
                             a::Float64, m::Float64)
    κ₁, κ₂ = 1.0 + Δκ, 1.0 - Δκ
    scale = 1 / sqrt(m)
    scale .* Dispersion.omega_op.(k_reduced .* a, κ₁, κ₂),
    scale .* Dispersion.omega_ac.(k_reduced .* a, κ₁, κ₂)
end

function plot_dispersion(Δκ_values::Vector{Float64},
                         k_reduced::AbstractVector,
                         a::Float64, m::Float64)
    freq_data = map(Δκ -> compute_frequencies(Δκ, k_reduced, a, m),
                    Δκ_values) # Mapping on

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
        ls  = cyc(LINESTYLES, i) 
        col = cyc(PALETTE_DELTA, i)
        lw  = cyc(LINEWIDTHS, i)

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