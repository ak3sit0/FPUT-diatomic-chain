module FPUTAnalysis

using LinearAlgebra, Statistics

export compute_modal_energies, sliding_window_avg, spectral_entropy

# ── Layer 1: Coordinate Transformations ──

"""
    compute_modal_energies(q, v, frequencies, V, m) 
Transforms physical coordinates to modal space and computes energy.
"""
function compute_modal_energies(q, v, freq, V, m)
    # U = M^(-1/2) * V
    U = Diagonal(1 ./ sqrt.(m)) * V
    # x_modal = U^-1 * q_physical
    q_m = U \ q
    v_m = U \ v
    return 0.5 .* (v_m.^2 .+ (freq.^2) .* q_m.^2)
end

# ── Layer 2: Signal Processing (Sliding Averages) ──

"""
    sliding_window_avg(series::AbstractVector, delta::Float64)
SICP-like: computes a moving average where window grows with time (Δ*t).
"""
function sliding_window_avg(E::AbstractVector{T}, Δ::Float64) where T
    n = length(E)
    S = pushfirst!(cumsum(E), zero(T))
    avgs = Vector{Float64}(undef, n)
    for t in 1:n
        start_t = max(1, Int(floor(Δ * t)))
        window = t - start_t + 1
        @inbounds avgs[t] = (S[t+1] - S[start_t]) / window
    end
    return avgs
end

# ── Layer 3: High-level Metrics ──

"""
    spectral_entropy(energies::Matrix, delta::Float64)
"""
function spectral_entropy(modal_E::Matrix, delta::Float64)
    # Apply sliding average to each mode (row)
    E_avg = stack(sliding_window_avg(row, delta) for row in eachrow(modal_E))'
    
    # Probabilities p_k = E_k / sum(E)
    total_E = sum(E_avg, dims=1)
    p = E_avg ./ (total_E .+ 1e-15)
    
    # S = -∑ p log p
    entropy = -vec(sum(p .* log.(p .+ 1e-15), dims=1))
    return entropy
end

end # module
