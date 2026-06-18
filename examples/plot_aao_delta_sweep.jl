using LinearAlgebra, CairoMakie, LaTeXStrings

include(joinpath(@__DIR__, "plot_gamma_with_resonance.jl"))

function plot_aao_delta_sweep(delta_values, alfa; Nk=601, Ngrid=201)
    kplot = collect(range(-π, π, length=Ngrid))
    aao_idx = 2  # "aao" es el índice 2 en names

    # Pre-calcular todos los |Γ_aao|
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
                  xlabelsize=16, ylabelsize=16, titlesize=16)

        last_hm = heatmap!(ax, kplot, kplot, Gaao;
                           colorrange=(0, global_max), colormap=:Blues)

        # Curva de resonancia aao: s1=1(a), s2=1(a), s3=2(o)
        D = resonance_matrix(kplot, kA, kB, 1, 1, 2)
        contour!(ax, kplot, kplot, D; levels=[0.0], color=:orangered, linewidth=2.5)
    end

    # Barra de color compartida
    Colorbar(fig[1:2, 4], last_hm;
             label=L"|\Gamma_{aao}(k_1,k_2)|", labelsize=16)

    # Supertítulo
    Label(fig[0, :],
          latexstring("\\left|\\Gamma_{aao}(k_1,k_2)\\right|,\\quad k_3=-k_1-k_2,\\quad \\alpha=$(alfa)"),
          fontsize=18)

    mkpath("results/figures/coupling")
    save("results/figures/coupling/Gamma_aao_delta_sweep.png", fig)
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_aao_delta_sweep([0.1, 0.2, 0.3, 0.4, 0.5], 0.1)
end
