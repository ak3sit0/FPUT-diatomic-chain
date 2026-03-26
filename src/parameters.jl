module Parameters

export DT, make_k_m

# ===========================================================================
# PHYSICAL CONSTANTS
# Values that are truly invariant across all experiments.
# Everything else (N, α, β, TMAX, …) belongs in a configs/*.toml file.
# ===========================================================================

"""Integration time step for the Störmer–Verlet scheme."""
const DT = 0.05

# ===========================================================================
# SYSTEM BUILDER
# ===========================================================================

"""
    make_k_m(N, DeltaK, DeltaM, boundary) → (k, m)

Constructs alternating spring-constant and mass vectors for the FPUT chain.

- `boundary = :fixed`    → length(k) = N+1
- `boundary = :periodic` → length(k) = N

Pattern: k[i] = 1 ± DeltaK,  m[i] = 1 ± DeltaM  (alternating sign).
"""
function make_k_m(N::Int, DeltaK::Float64, DeltaM::Float64, boundary::Symbol)
    if boundary == :fixed
        k = [1.0 + DeltaK * (-1.0)^(i-1) for i in 1:(N+1)]
    elseif boundary == :periodic
        k = [1.0 + DeltaK * (-1.0)^(i-1) for i in 1:N]
    else
        error("boundary must be :fixed or :periodic")
    end
    m = [1.0 + DeltaM * (-1.0)^(i-1) for i in 1:N]
    return k, m
end

end
