"""
    coupling_coefficients.jl

Rejilla 3×3 de |Γ_σ1σ2σ3(k₁,k₂)| para las 8 combinaciones de rama, con
k₃ = -k₁-k₂. Con `resonance=true` superpone además la curva de resonancia
propia de cada panel (antes `plot_gamma_with_resonance.jl`).

Usage:
  julia --project=. examples/coupling/coupling_coefficients.jl
"""

using LinearAlgebra, CairoMakie, LaTeXStrings
include("../../src/fput_coupling.jl");   using .FPUTCoupling
include("../../src/fput_analysis.jl");   using .FPUTAnalysis
include("../../src/plotting_utils.jl");  using .PlottingUtils

const BRANCH_IDX = Dict('a' => 1, 'o' => 2)

"""
    plot_gamma(kA, kB, alfa; resonance=false, Nk=601, Ngrid=201) -> Figure

Rejilla de los 8 |Γ|. `resonance=true` añade el contorno de resonancia por panel
y una barra de color por panel (variante PNG); `false` usa una única barra
compartida (variante PDF).
"""
function plot_gamma(kA, kB, alfa; resonance::Bool=false, Nk=601, Ngrid=201)
    kplot, names, Gamma = compute_gamma(kA, kB, alfa; Nk=Nk, Ngrid=Ngrid)
    maxabs = maximum(maximum(abs.(G)) for G in Gamma)

    # La variante con resonancia deja sitio a 8 barras de color, así que usa
    # tipografía algo menor que la de barra compartida.
    fs_label, fs_title, fs_super = resonance ? (22, 20, 22) : (26, 32, 34)

    fig = Figure(size=(1400, 1050))
    positions = [(i, j) for i in 1:3 for j in 1:3]

    local last_hm
    for n in 1:8
        row, col = positions[n]
        ax = Axis(fig[row, col],
                  xlabel=L"k_1", ylabel=L"k_2",
                  title=latexstring("\\left|\\Gamma_{$(names[n])}(k_1,k_2)\\right|"),
                  xlabelsize=fs_label, ylabelsize=fs_label, titlesize=fs_title)

        last_hm = heatmap!(ax, kplot, kplot, abs.(Gamma[n]);
                           colorrange=(0, maxabs), colormap=PALETTE_COUPLING,
                           rasterize=4)

        if resonance
            Colorbar(fig[row, col][1, 2], last_hm, width=30, labelsize=14)
            s1, s2, s3 = (BRANCH_IDX[names[n][1]], BRANCH_IDX[names[n][2]], BRANCH_IDX[names[n][3]])
            D = resonance_matrix(kplot, kA, kB, s1, s2, s3)
            contour!(ax, kplot, kplot, D; levels=[0.0], color=RESONANCE_LINE, linewidth=2.5)
        end
    end

    resonance || Colorbar(fig[1:3, 4], last_hm, width=32, labelsize=36)

    ax_empty = Axis(fig[3, 3]); hidedecorations!(ax_empty); hidespines!(ax_empty)

    Label(fig[0, :],
          latexstring("\\left|\\Gamma_{\\sigma_1\\sigma_2\\sigma_3}(k_1,k_2)\\right|,\\quad k_3=-k_1-k_2,\\quad \\kappa_A=$(kA),\\;\\kappa_B=$(kB)"),
          fontsize=fs_super)

    mkpath("results/figures/coupling")
    stem = resonance ? "Gamma_resonance_kA$(kA)_kB$(kB).png" : "Gamma_kA$(kA)_kB$(kB).pdf"
    save(joinpath("results/figures/coupling", stem), fig)
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_gamma(1.3, 0.7, 0.1)
    plot_gamma(1.1, 0.9, 0.1; resonance=true)
end
