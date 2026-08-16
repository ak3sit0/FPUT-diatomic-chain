module FPUTAnalysis

using LinearAlgebra, Statistics

export compute_modal_energies, sliding_window_avg, spectral_entropy

# ── Layer 1: Coordinate Transformations ──

"""
    compute_modal_energies(q, v, frequencies, V, m) 
Transforms physical coordinates to modal space and computes energy.
"""
function compute_modal_energies(q, v, freq, V, m)
    # Use mass-weighted coordinates explicitly for numerical clarity:
    # x = sqrt(m) .* q, v_x = sqrt(m) .* v
    x = sqrt.(m) .* q
    v_x = sqrt.(m) .* v
    # modal coordinates: project onto orthonormal eigenvectors V
    q_m = V' * x
    v_m = V' * v_x
    return 0.5 .* (v_m.^2 .+ (freq.^2) .* q_m.^2)
end

# ── Layer 2: Signal Processing (Sliding Averages) ──

"""
    sliding_window_avg!(out, S, E, delta)
In-place core: `S` is a caller-provided cumsum buffer of length `length(E)+1`,
reused across calls (e.g. one per mode in `spectral_entropy`) to avoid
reallocating it every time.
"""
function sliding_window_avg!(out::AbstractVector, S::AbstractVector, E::AbstractVector, Δ::Float64)
    n = length(E)
    @inbounds begin
        S[1] = zero(eltype(S))
        for i in 1:n
            S[i+1] = S[i] + E[i]
        end
        for t in 1:n
            start_t = max(1, Int(floor(Δ * t)))
            out[t] = (S[t+1] - S[start_t]) / (t - start_t + 1)
        end
    end
    return out
end

"""
    sliding_window_avg(series::AbstractVector, delta::Float64)
    computes a moving average where window grows with time (Δ*t).
"""
function sliding_window_avg(E::AbstractVector{T}, Δ::Float64) where T
    n = length(E)
    S = Vector{T}(undef, n + 1)
    avgs = Vector{Float64}(undef, n)
    sliding_window_avg!(avgs, S, E, Δ)
    return avgs
end

# ── Layer 3: High-level Metrics ──

"""
    spectral_entropy(energies::Matrix, delta::Float64)
"""
function spectral_entropy(modal_E::Matrix, delta::Float64)
    N, nt = size(modal_E)
    # Apply sliding average to each mode (row). E_avg is built directly in
    # (N, nt) orientation and S is reused across modes 
    E_avg = Matrix{Float64}(undef, N, nt)
    S = Vector{Float64}(undef, nt + 1)
    for i in 1:N
        sliding_window_avg!(view(E_avg, i, :), S, view(modal_E, i, :), delta)
    end

    # Probabilities p_k = E_k / sum(E)
    total_E = sum(E_avg, dims=1)
    p = E_avg ./ (total_E .+ 1e-15)
    
    # S = -∑ p log p
    entropy = -vec(sum(p .* log.(p .+ 1e-15), dims=1))
    return entropy
end

end # module
