using LinearAlgebra, CairoMakie, LaTeXStrings
include("../../src/fput_coupling.jl")
using .FPUTCoupling

function plot_gamma(kA, kB, alfa; Nk=601, Ngrid=201)
    kplot, names, Gamma = compute_gamma(kA, kB, alfa; Nk=Nk, Ngrid=Ngrid)

    maxabs = maximum(maximum(abs.(G)) for G in Gamma)

    fig = Figure(size=(1400, 1050))
    positions = [(i,j) for i in 1:3 for j in 1:3]

    local last_hm
    for n in 1:8
        row, col = positions[n]
        ax = Axis(fig[row, col],
                  xlabel=L"k_1", ylabel=L"k_2",
                  title=latexstring("\\left|\\Gamma_{$(names[n])}(k_1,k_2)\\right|"),
                  xlabelsize=26, ylabelsize=26, titlesize=32)
        last_hm = heatmap!(ax, kplot, kplot, abs.(Gamma[n]),
                      colorrange=(0, maxabs), colormap=:Blues,
                      rasterize=4)
    end

    # Barra de color compartida
    Colorbar(fig[1:3, 4], last_hm, width=32, labelsize=36)

    # empty tile
    ax_empty = Axis(fig[3, 3]); hidedecorations!(ax_empty); hidespines!(ax_empty)

    Label(fig[0, :],
          latexstring("\\left|\\Gamma_{\\sigma_1\\sigma_2\\sigma_3}(k_1,k_2)\\right|,\\quad k_3=-k_1-k_2,\\quad \\kappa_A=$(kA),\\;\\kappa_B=$(kB)"),
          fontsize=34)

    save("./results/figures/coupling/Gamma_kA$(kA)_kB$(kB).pdf", fig)
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_gamma(1.3, 0.7, 0.1)
end