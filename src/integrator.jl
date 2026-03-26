module Integrator

include("fput_equations.jl")
using .FPUTEquations
using DifferentialEquations

export run_simulation

function run_simulation(v0, q0, tspan, k, m, alpha, beta, DT, forces; save_every::Int=100)
    prob = SecondOrderODEProblem(FPUTEquations.fput_generalized, v0, q0, tspan, (k, m, alpha, beta, forces))
    # Save solution every save_every time steps to balance memory and resolution
    saveat = tspan[1]:(save_every*DT):tspan[2]
    return solve(prob, KahanLi8(), dt=DT, saveat=saveat)
end

end
