module FPUTCore

using LinearAlgebra

export SystemParams, make_system, fput_forces!, find_normal_modes

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
    fput_forces!(dv, q, p_ode, t)
The core Hamiltonian derivative. p_ode = (k, m, alpha, beta, boundary)
"""
function fput_forces!(dv, v, q, p_ode, t)
    k, m, alpha, beta, boundary = p_ode
    N = length(m)
    
    if boundary == :fixed
        # Handle interiors
        for i in 2:N-1
            dl, dr = q[i] - q[i-1], q[i+1] - q[i]
            dv[i] = (k[i+1]*dr - k[i]*dl + alpha*(k[i+1]*dr^2 - k[i]*dl^2) + beta*(k[i+1]*dr^3 - k[i]*dl^3)) / m[i]
        end
        # Boundaries
        dv[1] = (k[2]*(q[2]-q[1]) - k[1]*q[1] + alpha*(k[2]*(q[2]-q[1])^2 - k[1]*q[1]^2) + beta*(k[2]*(q[2]-q[1])^3 - k[1]*q[1]^3)) / m[1]
        
        dl_N, dr_N = q[N] - q[N-1], -q[N]
        dv[N] = (k[N+1]*dr_N - k[N]*dl_N + alpha*(k[N+1]*dr_N^2 - k[N]*dl_N^2) + beta*(k[N+1]*dr_N^3 - k[N]*dl_N^3)) / m[N]
        
    elseif boundary == :periodic
        for i in 1:N
            idx_prev = (i == 1) ? N : i - 1
            idx_next = (i == N) ? 1 : i + 1
            rl, rr = q[i] - q[idx_prev], q[idx_next] - q[i]

            # Forces derived from potential with k-scaling for nonlinear terms:
            # V_bond = k/2 * Δ^2 + α*k/3 * Δ^3 + β*k/4 * Δ^4
            # => F_bond = k*Δ + α*k*Δ^2 + β*k*Δ^3
            f_left  = k[idx_prev]*rl + alpha * k[idx_prev] * rl^2 + beta * k[idx_prev] * rl^3
            f_right = k[i]*rr + alpha * k[i] * rr^2 + beta * k[i] * rr^3
            dv[i] = (f_right - f_left) / m[i]
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
    λ, V = eigen(D)
    return sqrt.(max.(0.0, real.(λ))), V   # max evita sqrt de negativos por redondeo flotante
end

end # module FPUTCore
