module DynamicalMatrix

export find_normal_modes, build_dynamical_matrix

using LinearAlgebra

function build_dynamical_matrix(k::Vector{Float64}, m::Vector{Float64})
    N = length(m)
    D = zeros(N, N)

    if length(k) == N + 1 # Fixed boundary conditions
        for i in 1:N
            D[i,i] += (k[i]+ k[i+1]) / m[i]
            if i > 1
                D[i, i-1] = -k[i] / sqrt(m[i] * m[i-1])
            end
            if i < N
                D[i, i+1] = -k[i+1] / sqrt(m[i] * m[i+1])
            end
        end
    elseif length(k) == N # Periodic boundary conditions
        for i in 1:N
            D[i,i] += (k[i]+ k[mod1(i+1, N)]) / m[i]
            D[i, mod1(i-1, N)] = -k[i] / sqrt(m[i] * m[mod1(i-1, N)])
            D[i, mod1(i+1, N)] = -k[mod1(i+1, N)] / sqrt(m[i] * m[mod1(i+1, N)])
        end
    else
        error("Unsupported boundary conditions or inconsistent k and m lengths.")
    
    end

    return D
end

function find_normal_modes(k::Vector{Float64}, m::Vector{Float64})
    D = build_dynamical_matrix(k, m)
    λ, V = eigen(D)
    frequencies = sqrt.(abs.(λ))
    return frequencies, V
end

end
