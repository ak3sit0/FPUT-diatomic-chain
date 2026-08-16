using LinearAlgebra, CairoMakie, LaTeXStrings
include("../../src/fput_coupling.jl")
using .FPUTCoupling

# ── Plot principal ────────────────────────────────────────────────────────────

function plot_gamma_with_resonance(kA, kB, alfa; Nk=601, Ngrid=201)
    kplot, names, Gamma = compute_gamma(kA, kB, alfa; Nk=Nk, Ngrid=Ngrid)

    # Mapeo de nombre a índice de rama
    branch_idx = Dict('a' => 1, 'o' => 2)
    branch_by_panel = [(branch_idx[names[n][1]],
                        branch_idx[names[n][2]],
                        branch_idx[names[n][3]]) for n in 1:8]

    maxabs = maximum(maximum(abs.(G)) for G in Gamma)

    fig = Figure(size=(1400, 1050))
    positions = [(i, j) for i in 1:3 for j in 1:3]

    for n in 1:8
        row, col = positions[n]
        ax = Axis(fig[row, col],
                  xlabel=L"k_1", ylabel=L"k_2",
                  title=latexstring("\\left|\\Gamma_{$(names[n])}(k_1,k_2)\\right|"),
                  xlabelsize=22, ylabelsize=22, titlesize=20)

        # Heatmap del acoplamiento
        hm = heatmap!(ax, kplot, kplot, abs.(Gamma[n]),
                      colorrange=(0, maxabs), colormap=:Blues)
        Colorbar(fig[row, col][1, 2], hm, width=30, labelsize=14)

        # Curva de resonancia específica del panel
        s1, s2, s3 = branch_by_panel[n]
        D = resonance_matrix(kplot, kA, kB, s1, s2, s3)
        contour!(ax, kplot, kplot, D; levels=[0.0], color=:orangered, linewidth=2.5)
    end

    # Panel vacío (posición 3,3)
    ax_empty = Axis(fig[3, 3]); hidedecorations!(ax_empty); hidespines!(ax_empty)

    # Supertítulo
    Label(fig[0, :],
          latexstring("\\left\\|\\Gamma_{\\sigma_1\\sigma_2\\sigma_3}(k_1,k_2)\\right\\|,\\quad k_3=-k_1-k_2,\\quad \\kappa_A=$(kA),\\;\\kappa_B=$(kB)"),
          fontsize=22)

    mkpath("results/figures/coupling")
    save("./results/figures/coupling/Gamma_resonance_kA$(kA)_kB$(kB).png", fig)
    return fig
end

# ── Ejecución ─────────────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    plot_gamma_with_resonance(1.1, 0.9, 0.1)
end
