module FPUTCoupling

using LinearAlgebra

export wrap_pi, compute_eigenvectors, fetch_e, compute_gamma, omega_branch, resonance_matrix

# ── Utility: wrap k to [-π, π] ──
wrap_pi(k::Real) = mod(k + π, 2π) - π

# ── Eigenvectors and frequencies for diatomic chain ──
"""
    compute_eigenvectors(k_sample, kA, kB) -> (Estore, omega)

Compute eigenvectors and eigenfrequencies for diatomic chain with spring constants kA, kB.
Returns:
- Estore[ik, component, branch]: eigenvector components (2×2 for each k)
- omega[ik, branch]: frequencies (acoustic=1, optical=2)
"""
function compute_eigenvectors(k_sample::AbstractVector, kA::Float64, kB::Float64)
    Nk = length(k_sample)
    Estore = zeros(ComplexF64, Nk, 2, 2)
    omega  = zeros(Float64, Nk, 2)

    @inbounds for (ik, k) in enumerate(k_sample)
        S    = kA + kB * exp(+1im * k)
        shat = S / abs(S)
        Estore[ik, :, 1] = (1/√2) .* [-shat,  1.0]   # acoustic
        Estore[ik, :, 2] = (1/√2) .* [ shat,  1.0]   # optical

        ksum  = kA + kB
        kprod = kA * kB
        disc  = ksum^2 - 4kprod * sin(k/2)^2
        omega[ik, 1] = sqrt(max(ksum - sqrt(max(disc, 0.0)), 0.0))
        omega[ik, 2] = sqrt(ksum + sqrt(max(disc, 0.0)))
    end
    return Estore, omega
end

# ── Linear interpolation of eigenvectors on k-grid ──
"""
    fetch_e(k, k_sample, Estore) -> E

Linear interpolation of eigenvector E(k) from stored grid.
Returns 2×2 matrix with columns = [acoustic, optical] branches.
"""
function fetch_e(k::Real, k_sample::AbstractVector, Estore::Array{ComplexF64,3})
    k = wrap_pi(k)
    dk   = k_sample[2] - k_sample[1]
    idx  = clamp(searchsortedfirst(k_sample, k) - 1, 1, length(k_sample)-1)
    t    = (k - k_sample[idx]) / dk
    E    = zeros(ComplexF64, 2, 2)

    @inbounds for b in 1:2, c in 1:2
        E[c, b] = (1-t)*Estore[idx, c, b] + t*Estore[idx+1, c, b]
    end
    return E
end

# ── Coupling coefficients Γ_σ1σ2σ3 ──
"""
    compute_gamma(kA, kB, alfa; Nk=601, Ngrid=201) -> (kplot, names, Gamma)

Compute three-wave coupling coefficients |Γ_σ1σ2σ3(k1,k2)| on a 2D grid.
Returns:
- kplot: wavenumber grid for plotting
- names: branch labels (e.g., "aao", "ooa", etc.)
- Gamma: vector of 8 matrices, one per branch combination
"""
function compute_gamma(kA::Float64, kB::Float64, alfa::Float64; Nk::Int=601, Ngrid::Int=201)
    betaA = alfa * kA
    betaB = alfa * kB

    k_sample = collect(range(-π, π, length=Nk))
    Estore, _ = compute_eigenvectors(k_sample, kA, kB)

    kplot = collect(range(-π, π, length=Ngrid))

    names = ["aaa","aao","aoa","oaa","aoo","oao","ooa","ooo"]
    Gamma = [zeros(ComplexF64, Ngrid, Ngrid) for _ in 1:8]

    @inbounds for i1 in 1:Ngrid
        k1 = kplot[i1]
        for i2 in 1:Ngrid
            k2 = kplot[i2]
            k3 = wrap_pi(-k1 - k2)

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

# ── Dispersion relation (omega as function of k, branch) ──
"""
    omega_branch(k, kA, kB, branch::Int) -> ω

Acoustic (branch=1) or optical (branch=2) frequency at wavenumber k.
For diatomic chain with spring constants kA, kB.
"""
function omega_branch(k::Real, kA::Float64, kB::Float64, branch::Int)
    k = wrap_pi(k)
    ksum  = kA + kB
    kprod = kA * kB
    disc  = ksum^2 - 4kprod * sin(k/2)^2

    if branch == 1
        sqrt(max(ksum - sqrt(max(disc, 0.0)), 0.0))
    else
        sqrt(ksum + sqrt(max(disc, 0.0)))
    end
end

# ── Resonance matrix ──
"""
    resonance_matrix(kplot, kA, kB, s1, s2, s3) -> D

Residual of resonance condition: ω_s3(k3) - ω_s1(k1) - ω_s2(k2),
where k3 = -k1-k2. Returns Ngrid × Ngrid matrix with D[j,i] = residual(kplot[i], kplot[j]).
Note: D[j,i] corresponds to (k1=kplot[i], k2=kplot[j]).
"""
function resonance_matrix(kplot::AbstractVector, kA::Float64, kB::Float64, s1::Int, s2::Int, s3::Int)
    Ngrid = length(kplot)
    D = zeros(Ngrid, Ngrid)

    @inbounds for i in 1:Ngrid, j in 1:Ngrid
        k1 = kplot[i]
        k2 = kplot[j]
        k3 = wrap_pi(-k1 - k2)
        D[j, i] = omega_branch(k3, kA, kB, s3) -
                  omega_branch(k1, kA, kB, s1) -
                  omega_branch(k2, kA, kB, s2)
    end

    return D
end

end # module FPUTCoupling
