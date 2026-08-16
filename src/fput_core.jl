module FPUTCore

using LinearAlgebra

export SystemParams, make_system, fput_forces!, find_normal_modes, bond_potential

# ── Layer 1: Domain types (SICP: Data Abstraction) ──
struct SystemParams
    N::Int
    delta_k::Float64
    delta_m::Float64
    alpha::Float64
    beta::Float64
    boundary::Symbol # :fixed or :periodic
end

# ── Layer 2: System Construction ──
"""
    make_system(p::SystemParams) -> (k, m)
Constructs the chain's physical properties. Pattern: 1 ± Δ (alternating).
"""
function make_system(p::SystemParams)
    k_len = p.boundary == :fixed ? p.N + 1 : p.N
    k = [1.0 + p.delta_k * (-1.0)^(i-1) for i in 1:k_len]
    m = [1.0 + p.delta_m * (-1.0)^(i-1) for i in 1:p.N]
    return k, m
end
# ── Layer 3: Physics Kernels (Functional style within mutation) ──
"""
    bond_force(kj, d, alpha, beta) -> F

Bond force j with elongation `d`. From potential with k scaling:
  V_bond = k/2·Δ² + α·k/3·Δ³ + β·k/4·Δ⁴   ⇒   F_bond = k·Δ + α·k·Δ² + β·k·Δ³
Factored as k·Δ·(1 + Δ·(α + β·Δ)) to save multiplications.
"""
@inline bond_force(kj, d, alpha, beta) = kj * d * (1 + d * (alpha + beta * d)) # @inline for performance, since this is called in tight loops.
@inline bond_force_a(kj, d, alpha)     = kj * d * (1 + alpha * d)   # case β = 0

"""
    bond_potential(kj, d, alpha, beta) -> V

Bond potential energy: V_bond = k/2·Δ² + α·k/3·Δ³ + β·k/4·Δ⁴
Used for energy conservation checks in tests.
"""
@inline bond_potential(kj, d, alpha, beta) = kj * (0.5*d^2 + (alpha/3)*d^3 + (beta/4)*d^4)

# Forces by BOND (not by site): each bond is evaluated only once, then
# dv[i] = (F_right - F_left)·inv_m[i]. Previous version computed each bond
# twice (as f_right of i and as f_left of i+1) and had boundary conditionals
# inside the loop, preventing vectorization.

function _bonds_periodic!(F, q, k, alpha, beta, N)
    @inbounds if beta == 0.0
        @simd for i in 1:N-1
            F[i] = bond_force_a(k[i], q[i+1] - q[i], alpha)
        end
        F[N] = bond_force_a(k[N], q[1] - q[N], alpha)
    else
        @simd for i in 1:N-1
            F[i] = bond_force(k[i], q[i+1] - q[i], alpha, beta)
        end
        F[N] = bond_force(k[N], q[1] - q[N], alpha, beta)
    end
end

function _bonds_fixed!(F, q, k, alpha, beta, N)
    # Bond 1: wall–site 1 (d = q₁);  bond N+1: site N–wall (d = −q_N)
    @inbounds if beta == 0.0
        F[1] = bond_force_a(k[1], q[1], alpha)
        @simd for j in 2:N
            F[j] = bond_force_a(k[j], q[j] - q[j-1], alpha)
        end
        F[N+1] = bond_force_a(k[N+1], -q[N], alpha)
    else
        F[1] = bond_force(k[1], q[1], alpha, beta)
        @simd for j in 2:N
            F[j] = bond_force(k[j], q[j] - q[j-1], alpha, beta)
        end
        F[N+1] = bond_force(k[N+1], -q[N], alpha, beta)
    end
end

"""
    fput_forces!(dv, v, q, p_ode, t)
The core Hamiltonian derivative. p_ode = (k, inv_m, alpha, beta, boundary, F)
`inv_m` = 1 ./ m precomputed; `F` = force buffer per bond (length N for
:periodic, N+1 for :fixed). The buffer is per-call to `solve_fput`, so each
thread has its own.
"""
function fput_forces!(dv, v, q, p_ode, t)
    k, inv_m, alpha, beta, boundary, F = p_ode
    N = length(inv_m)

    if boundary == :periodic
        _bonds_periodic!(F, q, k, alpha, beta, N)
        @inbounds begin
            dv[1] = (F[1] - F[N]) * inv_m[1]
            @simd for i in 2:N
                dv[i] = (F[i] - F[i-1]) * inv_m[i]
            end
        end
    elseif boundary == :fixed
        _bonds_fixed!(F, q, k, alpha, beta, N)
        @inbounds @simd for i in 1:N
            dv[i] = (F[i+1] - F[i]) * inv_m[i]
        end
    end
end

# ── Layer 4: Spectral Theory (Normal Modes) ──

"""
    find_normal_modes(k, m, boundary) -> (frequencies, V)
V is the matrix of eigenvectors.
"""
function find_normal_modes(k, m, boundary)
    N = length(m)
    D = zeros(N, N)
    if boundary == :fixed
        for i in 1:N
            D[i,i] += (k[i]+ k[i+1]) / m[i]
            i > 1 && (D[i, i-1] = -k[i] / sqrt(m[i]*m[i-1]))
            i < N && (D[i, i+1] = -k[i+1] / sqrt(m[i]*m[i+1]))
        end
    else # Periodic
        for i in 1:N
            prev = mod1(i-1, N)
            nxt  = mod1(i+1, N)
            # diagonal: sum of left and right spring constants divided by mass
            D[i,i] += (k[prev] + k[i]) / m[i]
            # coupling to previous mass via left spring k[prev]
            D[i, prev] = -k[prev] / sqrt(m[i] * m[prev])
            # coupling to next mass via right spring k[i]
            D[i, nxt] = -k[i] / sqrt(m[i] * m[nxt])
        end
    end
    # Symmetric(D) selects LAPACK's symmetric eigensolver (syevr): faster than the
    # general dense eigen(), and guarantees real eigenvalues (no real.() needed).
    λ, V = eigen(Symmetric(D))
    return sqrt.(max.(0.0, λ)), V   # max avoids sqrt of tiny negative λ from rounding
end

end # module FPUTCore
