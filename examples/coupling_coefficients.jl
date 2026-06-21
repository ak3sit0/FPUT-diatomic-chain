using LinearAlgebra, CairoMakie, LaTeXStrings

function compute_eigenvectors(k_sample, kA, kB)
    Nk = length(k_sample)
    # Estore[ik, component, branch]: branch 1=acoustic, 2=optical
    Estore = zeros(ComplexF64, Nk, 2, 2)
    omega  = zeros(Float64, Nk, 2)
    for (ik, k) in enumerate(k_sample)
        S    = kA + kB * exp(+1im * k)
        shat = S / abs(S)
        Estore[ik, :, 1] = (1/√2) .* [-shat,  1.0]   # acoustic
        Estore[ik, :, 2] = (1/√2) .* [ shat,  1.0]   # optical
        ksum  = kA + kB
        kprod = kA * kB
        disc  = ksum^2 - 4kprod * sin(k/2)^2
        omega[ik, 1] = sqrt(max((ksum - sqrt(disc)), 0.0))
        omega[ik, 2] = sqrt(ksum + sqrt(disc))
    end
    return Estore, omega
end

function fetch_e(k, k_sample, Estore)
    k = mod(k + π, 2π) - π
    # linear interpolation for each component and branch
    dk   = k_sample[2] - k_sample[1]
    idx  = clamp(searchsortedfirst(k_sample, k) - 1, 1, length(k_sample)-1)
    t    = (k - k_sample[idx]) / dk
    E    = zeros(ComplexF64, 2, 2)
    for b in 1:2, c in 1:2
        E[c, b] = (1-t)*Estore[idx, c, b] + t*Estore[idx+1, c, b]
    end
    return E  # E[:, 1]=acoustic, E[:, 2]=optical
end

function compute_gamma(kA, kB, alfa; Nk=601, Ngrid=201)
    betaA = alfa * kA
    betaB = alfa * kB

    k_sample = collect(range(-π, π, length=Nk))
    Estore, _ = compute_eigenvectors(k_sample, kA, kB)

    wrap = k -> mod(k + π, 2π) - π
    kplot = collect(range(-π, π, length=Ngrid))

    # branch index: 1=acoustic(a), 2=optical(o)
    # nameIdx = (s1-1)*4 + (s2-1)*2 + (s3-1) + 1
    names = ["aaa","aao","aoa","oaa","aoo","oao","ooa","ooo"]
    Gamma = [zeros(ComplexF64, Ngrid, Ngrid) for _ in 1:8]

    for i1 in 1:Ngrid
        k1 = kplot[i1]
        for i2 in 1:Ngrid
            k2 = kplot[i2]
            k3 = wrap(-k1 - k2)

            E1 = fetch_e(k1, k_sample, Estore)
            E2 = fetch_e(k2, k_sample, Estore)
            E3 = fetch_e(k3, k_sample, Estore)

            for s1 in 1:2, s2 in 1:2, s3 in 1:2
                e1 = E1[:, s1]
                e2 = E2[:, s2]
                e3 = E3[:, s3]

                DA1 = e1[2] - e1[1]
                DA2 = e2[2] - e2[1]
                DA3 = e3[2] - e3[1]

                DB1 = exp(+1im*k1)*e1[1] - e1[2]
                DB2 = exp(+1im*k2)*e2[1] - e2[2]
                DB3 = exp(+1im*k3)*e3[1] - e3[2]

                val = betaA*(DA1*DA2*DA3) + betaB*(DB1*DB2*DB3)
                idx = (s1-1)*4 + (s2-1)*2 + (s3-1) + 1
                Gamma[idx][i2, i1] = val
            end
        end
    end
    return kplot, names, Gamma
end

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
                  xlabelsize=22, ylabelsize=22, titlesize=28)
        last_hm = heatmap!(ax, kplot, kplot, abs.(Gamma[n]),
                      colorrange=(0, maxabs), colormap=:Blues)
    end

    # Barra de color compartida
    Colorbar(fig[1:3, 4], last_hm, width=32, labelsize=32)

    # empty tile
    ax_empty = Axis(fig[3, 3]); hidedecorations!(ax_empty); hidespines!(ax_empty)

    Label(fig[0, :],
          latexstring("\\left|\\Gamma_{\\sigma_1\\sigma_2\\sigma_3}(k_1,k_2)\\right|,\\quad k_3=-k_1-k_2,\\quad \\kappa_A=$(kA),\\;\\kappa_B=$(kB)"),
          fontsize=30)

    save("./results/figures/coupling/Gamma_kA$(kA)_kB$(kB).png", fig)
    return fig
end

# Ejecución
plot_gamma(1.2, 0.8, 0.2)