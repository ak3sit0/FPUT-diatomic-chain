using LinearAlgebra, CairoMakie, LaTeXStrings
include("../../src/fput_coupling.jl");   using .FPUTCoupling
include("../../src/fput_analysis.jl");   using .FPUTAnalysis
include("../../src/plotting_utils.jl");  using .PlottingUtils

function plot_aao_delta_sweep(delta_values, alfa; Nk=601, Ngrid=201)
    kplot = collect(range(-π, π, length=Ngrid)) # Range of values of k to plot
    aao_idx = 2  # "--+" is the second index

    # Precompute all |Γ_aao|
    aao_data = map(delta_values) do delta
        kA = 1.0 + delta
        kB = 1.0 - delta
        _, _, Gamma = compute_gamma(kA, kB, alfa; Nk=Nk, Ngrid=Ngrid)
        abs.(Gamma[aao_idx])
    end

    global_max = maximum(maximum, aao_data)

    fig = Figure(size=(1300, 800))

    local last_hm
    for (idx, (delta, Gaao)) in enumerate(zip(delta_values, aao_data))
        row = div(idx - 1, 3) + 1
        col = mod(idx - 1, 3) + 1
        kA  = 1.0 + delta
        kB  = 1.0 - delta

        ax = Axis(fig[row, col],
                  xlabel=L"k_1", ylabel=L"k_2",
                  title=latexstring("\\Delta\\kappa = $(delta)"),
                  xlabelsize=27, ylabelsize=27, titlesize=30)

        last_hm = heatmap!(ax, kplot, kplot, Gaao;
                           colorrange=(0, global_max), colormap=PALETTE_COUPLING,
                           rasterize=4)

        # aao resonance curve: s1=1(a), s2=1(a), s3=2(o)
        D = resonance_matrix(kplot, kA, kB, 1, 1, 2)
        contour!(ax, kplot, kplot, D; levels=[0.0], color=RESONANCE_LINE, linewidth=2.5)
    end

    # Shared colorbar
    Colorbar(fig[1:2, 4], last_hm;
             label=L"|\Gamma_{--+}(k_1,k_2)|", width=35, labelsize=40)

    # Supertitle
    Label(fig[0, :],
          latexstring("\\left|\\Gamma_{--+}(k_1,k_2)\\right|,\\quad k_3=-k_1-k_2,\\quad \\alpha=$(alfa)"),
          fontsize=34)

    mkpath("results/figures/coupling")
    save("results/figures/coupling/Gamma_aao_delta_sweep.pdf", fig)
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_aao_delta_sweep([0.05, 0.1, 0.2, 0.3, 0.4, 0.5], 0.1)
end
