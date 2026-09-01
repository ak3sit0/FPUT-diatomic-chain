module FPUTCore

using LinearAlgebra, Random

export SystemParams, make_system, fput_forces!, find_normal_modes, bond_potential,
       derive_band_indices, band_phase_ic, ref_frequency, OMEGA_MIN

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

# ── Layer 5: Band-selective initial conditions ──

"""
    derive_band_indices(branch, N, freq, boundary; k_band_start, k_band_end)

Returns Vector{Int} (1-based) of mode indices for the requested branch.
`freq` must be sorted ascending before calling this function.
If k_band_start/k_band_end are given, uses them directly (manual override).
"""
function derive_band_indices(branch::String, N::Int, freq::Vector{Float64},
                              boundary::Symbol;
                              k_band_start::Union{Int,Nothing}=nothing,
                              k_band_end::Union{Int,Nothing}=nothing)
    if !isnothing(k_band_start) && !isnothing(k_band_end)
        # The automatic branch below deliberately avoids mode 1 under PBC (k=0
        # translation, ω≈3e-8 from roundoff). When indices are given by hand that
        # protection did not apply: it is replicated here, plus the range bounds.
        1 <= k_band_start || error("k_band_start=$k_band_start must be ≥ 1 (1-based indices)")
        k_band_end <= N || error("k_band_end=$k_band_end exceeds N=$N")
        k_band_start <= k_band_end ||
            error("Empty band: k_band_start=$k_band_start > k_band_end=$k_band_end")
        boundary == :periodic && k_band_start == 1 && error(
            "k_band_start=1 with PBC includes the uniform translation k=0 (ω≈0, not a " *
            "physical mode): it would collapse scaled_t and make TMAX diverge. Use k_band_start ≥ 2.")
        return collect(k_band_start:k_band_end)
    end

    if boundary == :periodic
        isodd(N) && error("The diatomic chain with PBC requires even N (alternating springs); N=$N")

        # The acoustic/optical gap sits STRUCTURALLY at N/2: N/2 modes per branch.
        # It is NOT detected via argmax(diff(freq)): for small N and small Δκ the
        # intraband spacing (~1/N) exceeds the gap (~Δκ) and argmax lands inside the
        # acoustic branch. Measured: N=32 Δκ=0.10 gave gap_idx=3 (band 2:3 instead of
        # 2:16); N=64 Δκ=0.05 also failed. From N≥128 argmax gets it right, but there
        # is no need to rely on it when there is an exact answer.
        gap_idx = N ÷ 2

        # argmax is kept only as a diagnostic.
        argmax_idx = argmax(diff(freq)[2:end]) + 1
        if argmax_idx != gap_idx
            @warn "Gap by argmax ($argmax_idx) ≠ structural ($gap_idx): the intraband " *
                  "spacing exceeds the gap. Using the structural one (N/2)." N branch
        end

        branch == "acoustic" && return collect(2:gap_idx)   # excludes translation (ω≈0)
        branch == "optical"  && return collect(gap_idx+1:N)
        error("branch must be 'acoustic' or 'optical', got: '$branch'")
    else  # :fixed
        branch == "acoustic" && return collect(1:N)
        error("Fixed boundary has no distinguishable optical branch")
    end
end

"""
    band_phase_ic(k_band, E_total, N, freq, V, m, seed) -> (q0, v0)

Band-selective initial condition with uniform random phases.

E_total distributed uniformly over the modes in k_band:
  E_j = E_total / length(k_band)  for j ∈ k_band
  Q_j =  sqrt(2·E_j) / ω_j · cos(φ_j),  P_j = -sqrt(2·E_j) · sin(φ_j)
  q = (1/√m) · V · Q,   v = (1/√m) · V · P

Post-construction validation with @assert (machine precision).
"""
function band_phase_ic(k_band::AbstractVector{Int}, E_total::Float64, N::Int,
                       freq::Vector{Float64}, V::Matrix{Float64},
                       m::Vector{Float64}, seed::Int)
    rng    = Random.Xoshiro(seed)
    E_per  = E_total / length(k_band)
    Q      = zeros(N)
    P      = zeros(N)

    for j in k_band
        freq[j] < 1e-10 && continue   # Goldstone mode (translation, ω≈0)
        φ    = rand(rng) * 2π
        A    = sqrt(2 * E_per)
        Q[j] =  A / freq[j] * cos(φ)
        P[j] = -A            * sin(φ)
    end

    inv_sqrt_m = 1.0 ./ sqrt.(m)
    q0 = (V * Q) .* inv_sqrt_m
    v0 = (V * P) .* inv_sqrt_m

    # Project back to verify
    x        = sqrt.(m) .* q0
    vx       = sqrt.(m) .* v0
    Q_check  = V' * x
    P_check  = V' * vx
    E_check  = 0.5 .* (P_check.^2 .+ (freq.^2) .* Q_check.^2)
    E_out    = sum(E_check[setdiff(1:N, k_band)])
    E_in     = sum(E_check[k_band])
    tol      = 1e-8 * E_total

    @assert E_out < tol        "Energy leak outside the band: E_out=$(E_out) (tol=$(tol))"
    @assert abs(E_in - E_total) < tol "Incorrect normalization: E_in=$(E_in) vs E_total=$(E_total)"

    return q0, v0
end

# ── Layer 6: Reference frequency (time scale) ──

"""Below this, a frequency is numerical noise rather than a physical mode."""
const OMEGA_MIN = 1e-6

"""
    ref_frequency(freq, ref_mode) -> Float64

Frequency that sets the time scale (`scaled_t = t·ω/2π`) and, when a config gives
`scaled_t_max`, the integration budget (`TMAX = scaled_t_max·2π/ω`).

Errors when `ω ≈ 0`. Under PBC the mode 1 is the uniform k=0 translation, whose
eigenvalue is roundoff (`λ ~ 1e-16`, so `ω ~ 3e-8` — finite, not `Inf`), which makes
the failure silent: `scaled_t` collapses to ~0 and `TMAX` blows up by ~10⁶.
`ω₂ = 6.25/N`, so the threshold only bites for N ≳ 10⁶ — there is no physical mode
between it and the noise floor.
"""
function ref_frequency(freq::AbstractVector, ref_mode::Integer)
    1 <= ref_mode <= length(freq) ||
        error("ref_mode=$ref_mode out of range (1..$(length(freq)))")
    ω = freq[ref_mode]
    ω >= OMEGA_MIN || error(
        "ref_mode=$ref_mode has ω=$ω ≈ 0: it is the uniform k=0 translation of PBC, not a " *
        "physical mode. scaled_t would collapse and TMAX would diverge ~1e6×. With " *
        "boundary=:periodic use a mode ≥ 2.")
    ω
end

end # module FPUTCore
