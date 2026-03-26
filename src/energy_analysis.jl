module EnergyAnalysis

using LinearAlgebra, Plots, Statistics, NumericalIntegration

export find_energies, 
       kinetic_energy, potential_energy, total_energy, compute_energy,
       compute_entropy, 
       compute_eta, average_eta, cumulative_average

"""
    find_energies(v, q, frequencies, normal_matrix, m) -> N x T matrix

Computes modal energies for all time steps using the normalized coordinate transform.
"""
function find_energies(v, q, frequencies, normal_matrix, m)
    # Vectorized: transform todas las columnas de golpe
    sqrt_inv_mass = Diagonal(sqrt.(1 ./ m))
    U_matrix = sqrt_inv_mass * normal_matrix
    normal_q = U_matrix \ q         # N x T
    normal_v = U_matrix \ v         # N x T
    energies = 0.5 .* (normal_v .^ 2 .+ (frequencies .^ 2) .* (normal_q .^ 2))
    return energies
end

function kinetic_energy(v, m)
    return 0.5 * sum(m .* v.^2)
end

function potential_energy(q::Vector{Float64}, k::Vector{Float64})
    N = length(q)

    V = 0.0
    diff = q[2:end] .- q[1:end-1]


    if length(k) == N + 1 #Fixed boundary condition
        V = 0.5 * sum(k[2:end-1] .* diff.^2) # 
        V += 0.5 * (k[1] * q[1].^2 + k[end] * q[end].^2) # Boundary terms
    elseif length(k) == N #Periodic boundary condition 
        V += 0.5 * sum(k[1:N-1] .* diff.^2)
        V += 0.5 * k[N] * (q[N] - q[1])^2  # Wrap-around spring
    else
        error("Boundary condition not allowed")
    end

    return V
end

function total_energy(q, v, m, k)
    T = kinetic_energy(v, m)
    V = potential_energy(q, k)
    return T + V
end

function compute_energy(v, q, m, k)
    total_time = size(q, 2)
    energies = zeros(total_time)
    for t in 1:total_time
        energies[t] = total_energy(q[:, t], v[:, t], m, k)
    end
    return energies
end

function compute_entropy(energies::Matrix, delta::Float64=0.6)
    # Compute spectral entropy using sliding-delta averaging
    # delta controls window size: window_size = floor(delta * t)
    results = sliding_delta_average.(eachrow(energies), delta)
    averaged_energies = reduce(hcat, results)'

    # Probability distribution over modes
    p = averaged_energies ./ sum(averaged_energies, dims=1)
    p[p .< 1e-12] .= 1e-12 # Avoid log(0)
    S = -sum(p .* log.(p), dims=1)
    return S[1,:]
end

function compute_eta(spectral_entropy::Vector{Float64})
    S_max = maximum(spectral_entropy)
    S_0 = spectral_entropy[1] # Use the initial entropy, S(t=0)
    
    denominator = S_0 - S_max
    
    # If the entropy barely changes, the system is stable. Return zeros.
    if abs(denominator) < 1e-9
        return zeros(length(spectral_entropy))
    end
    
    # This is the correct formula from the MATLAB script
    return (spectral_entropy .- S_max) ./ denominator
end

function average_eta(eta)
    return cumsum(eta) ./ (1:length(eta))
end

"""
Calculates the cumulative average of a 2D matrix over its second dimension (time).
"""

function cumulative_average(series::Matrix)
    num_points = size(series, 2)
    avg_series = similar(series, Float64)
    current_sum = zeros(size(series, 1))
    for t in 1:num_points
        current_sum .+= series[:, t]
        avg_series[:, t] = current_sum ./ t
    end
    return avg_series
end

"""
Calculates the cumulative average of a 1D vector.
"""
function cumulative_average(series::Vector)
    num_points = length(series)
    avg_series = similar(series, Float64)
    current_sum = 0.0
    for t in 1:num_points
        current_sum += series[t]
        avg_series[t] = current_sum / t
    end
    return avg_series
end

"""
    windowed_average(series, window_size) -> averaged matrix

Computes moving average over window_size samples within each time slice.
Uses prefix sums for O(n) performance.
"""
function windowed_average(series::AbstractMatrix{<:Real}, window_size::Integer)
    @assert window_size >= 1 "window_size must be >= 1"
    m, n = size(series)
    half = div(window_size, 2)

    # prefix sums across time (columns)
    prefix = cumsum(Float64.(series), dims=2)  # Float64 for numeric stability
    out = similar(prefix)

    for j in 1:n
        a = max(1, j - half)
        b = min(n, j + half)
        len = b - a + 1
        if a == 1
            out[:, j] = prefix[:, b] ./ len
        else
            out[:, j] = (prefix[:, b] .- prefix[:, a-1]) ./ len
        end
    end

    return out
end


"""
    sliding_delta_average(E, Δ) -> averaging vector

Sliding-window average with window_size = floor(Δ * t) at time t.
Uses cumsum for O(n) performance.
"""
function sliding_delta_average(E::AbstractVector{T}, Δ::Float64) where T
    n = length(E)
    S = pushfirst!(cumsum(E), zero(T))
    avgs = Vector{Float64}(undef, n)
    
    for t in 1:n
        start_t = max(1, Int(floor(Δ * t)))
        current_window_width = t - start_t + 1
        @inbounds avgs[t] = (S[t+1] - S[start_t]) / current_window_width
    end
    return avgs
end

end

