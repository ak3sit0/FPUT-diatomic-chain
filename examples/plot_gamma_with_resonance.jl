using LinearAlgebra, CairoMakie, LaTeXStrings

# ── Eigenvectores y coeficientes Γ ─────────────────────────────────────────────

function compute_eigenvectors(k_sample, kA, kB)
    Nk = length(k_sample)
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
    dk   = k_sample[2] - k_sample[1]
    idx  = clamp(searchsortedfirst(k_sample, k) - 1, 1, length(k_sample)-1)
    t    = (k - k_sample[idx]) / dk
    E    = zeros(ComplexF64, 2, 2)
    for b in 1:2, c in 1:2
        E[c, b] = (1-t)*Estore[idx, c, b] + t*Estore[idx+1, c, b]
    end
    return E
end

function compute_gamma(kA, kB, alfa; Nk=601, Ngrid=201)
    betaA = alfa * kA
    betaB = alfa * kB

    k_sample = collect(range(-π, π, length=Nk))
    Estore, _ = compute_eigenvectors(k_sample, kA, kB)

    wrap = k -> mod(k + π, 2π) - π
    kplot = collect(range(-π, π, length=Ngrid))

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

# ── Dispersión para curvas de resonancia ──────────────────────────────────────

function omega_branch(k, kA, kB, branch)
    k = mod(k + π, 2π) - π
    ksum  = kA + kB
    kprod = kA * kB
    disc  = ksum^2 - 4kprod * sin(k/2)^2
    if branch == 1
        sqrt(max(ksum - sqrt(max(disc, 0.0)), 0.0))
    else
        sqrt(ksum + sqrt(max(disc, 0.0)))
    end
end

# ── Matriz residual de resonancia ─────────────────────────────────────────────

function resonance_matrix(kplot, kA, kB, s1, s2, s3)
    Ngrid = length(kplot)
    D = zeros(Ngrid, Ngrid)
    for i in 1:Ngrid, j in 1:Ngrid
        k1 = kplot[i]
        k2 = kplot[j]
        k3 = mod(-k1 - k2 + π, 2π) - π
        D[j,i] = omega_branch(k3, kA, kB, s3) -
                 omega_branch(k1, kA, kB, s1) -
                 omega_branch(k2, kA, kB, s2)
    end
    return D
end

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
        Colorbar(fig[row, col][1, 2], hm)

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
