module SimulationRunner

# Include all necessary modules
include("parameters.jl"); using .Parameters
include("dynamical_matrix.jl"); using .DynamicalMatrix
include("fput_equations.jl"); using .FPUTEquations
include("integrator.jl"); using .Integrator
include("energy_analysis.jl"); using .EnergyAnalysis

using LinearAlgebra

export run_and_analyze

# make_k_m lives in Parameters — imported above via `using .Parameters`

"""
    run_and_analyze(; N, TMAX, DT, alpha, beta, DeltaK, DeltaM, boundary, init_mode)

Ejecuta una simulación y devuelve (energies_per_mode, scaled_time, frequencies, v_final).
Nota: es pura en cuanto a I/O (ideal para threading).
"""
function run_and_analyze(; N::Int,
                          TMAX::Float64,
                          DT::Float64,
                          alpha::Float64=0.0,
                          beta::Float64=0.0,
                          DeltaK::Float64=0.0,
                          DeltaM::Float64=0.0,
                          boundary::Symbol=:fixed,
                          init_mode::Int=2,
                          energy::Float64=0.45,
                          q0_override::Union{Nothing,Vector{Float64}}=nothing,
                          v0_override::Union{Nothing,Vector{Float64}}=nothing,
                          tspan_override::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                          save_every::Int=100)
    # Build system
    k, m = make_k_m(N, DeltaK, DeltaM, boundary)

    # Normal modes / initial conditions
    frequencies, normal_matrix = find_normal_modes(k, m)
    
    # Use override state if provided, otherwise construct from init_mode
    if q0_override !== nothing && v0_override !== nothing
        q0 = q0_override
        v0 = v0_override
    else
        U = Diagonal(sqrt.(1 ./ m)) * normal_matrix
        # Calculate amplitude to satisfy the target energy
        # E = 0.5 * omega^2 * A^2  => A = sqrt(2E) / omega
        amplitude = sqrt(2 * energy) / frequencies[init_mode]
        q0 = amplitude * U[:, init_mode]
        v0 = zeros(N)
    end
    
    forces = zeros(N)

    # Time span: use override if provided, else (0.0, TMAX)
    tspan = tspan_override !== nothing ? tspan_override : (0.0, TMAX)

    # Run integrator with save_every
    sol = Integrator.run_simulation(v0, q0, tspan, k, m, alpha, beta, DT, forces; save_every=save_every)

    # Extract q, v efficiently
    Nt = length(sol.u)
    q = zeros(N, Nt)
    v = zeros(N, Nt)
    for (idx, u) in enumerate(sol.u)
        v[:, idx] = u.x[1]
        q[:, idx] = u.x[2]
    end

    # Modal energies and scaled time
    energies_per_mode = find_energies(v, q, frequencies, normal_matrix, m)
    # Scale time using angular frequency omega (rad/s): cycles = omega / (2π)
    scaled_time = sol.t .* frequencies[init_mode] ./ (2π)

    # Final state for checkpointing / continuation
    v_final = sol.u[end].x[1]
    q_final = sol.u[end].x[2]

    return energies_per_mode, scaled_time, frequencies, v_final, q_final
end

end
