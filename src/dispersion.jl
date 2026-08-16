module Dispersion

export omega_ac, omega_op, omega_compact_minus, omega_compact_plus, domega_dk

# ── Universal wrap to [-π, π] ──
wrap_pi(k::Real) = mod(k + π, 2π) - π

# ── Dispersion: absolute form with spring constants ──
"""
    omega_ac(k, kA, kB) -> ω

Acoustic branch frequency for diatomic chain with springs kA, kB.
"""
function omega_ac(k::Real, kA::Float64, kB::Float64)
    k = wrap_pi(k)
    ksum  = kA + kB
    kprod = kA * kB
    disc  = ksum^2 - 4kprod * sin(k/2)^2
    sqrt(max(ksum - sqrt(max(disc, 0.0)), 0.0))
end

"""
    omega_op(k, kA, kB) -> ω

Optical branch frequency for diatomic chain with springs kA, kB.
"""
function omega_op(k::Real, kA::Float64, kB::Float64)
    k = wrap_pi(k)
    ksum  = kA + kB
    kprod = kA * kB
    disc  = ksum^2 - 4kprod * sin(k/2)^2
    sqrt(ksum + sqrt(max(disc, 0.0)))
end

# ── Dispersion: compact form with normalized delta (κ* = 1 - Δκ²) ──
"""
    omega_compact_minus(k, delta) -> ω

Acoustic branch in compact parametrization: κ* = 1 - Δκ².
Valid for any convention A/B (κ₁ = 1+Δκ, κ₂ = 1-Δκ is equivalent).
"""
function omega_compact_minus(k::Real, delta::Float64)
    k = wrap_pi(k)
    A = sqrt(max(1 - (1 - delta^2) * sin(k/2)^2, 0.0))
    sqrt(max(2 - 2A, 0.0))
end

"""
    omega_compact_plus(k, delta) -> ω

Optical branch in compact parametrization: κ* = 1 - Δκ².
"""
function omega_compact_plus(k::Real, delta::Float64)
    k = wrap_pi(k)
    A = sqrt(max(1 - (1 - delta^2) * sin(k/2)^2, 0.0))
    sqrt(2 + 2A)
end

# ── Dispersion derivatives ──
"""
    domega_dk(k, kA, kB, branch) -> dω/dk

Analytic group velocity for acoustic (branch=1) or optical (branch=2) branch.
Derived from d(ω²)/dk = ∓ 2kA·kB·sin(k) / (2√disc).
"""
function domega_dk(k::Real, kA::Float64, kB::Float64, branch::Int)
    k = wrap_pi(k)
    ksum = kA + kB
    kprod = kA * kB
    sin2  = sin(k/2)^2
    disc  = ksum^2 - 4kprod * sin2
    disc  = max(disc, 1e-14)

    # d(disc)/dk = −2kA·kB·sin(k)
    d_disc_dk = -2kprod * sin(k)

    ω = (branch == 1) ? sqrt(max(ksum - sqrt(disc), 0.0)) :
                        sqrt(ksum + sqrt(disc))
    ω < 1e-12 && return 0.0

    sign = (branch == 1) ? -1.0 : +1.0
    return sign * d_disc_dk / (4ω * sqrt(disc))
end

end # module Dispersion
