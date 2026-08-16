module FPUTFastRunner

using DifferentialEquations
using ..FPUTCore # FPUTCore is the module that contains the core physics and system definitions

export solve_fput # Functional wrapper for the ODE solver.

"""
    solve_fput(p::SystemParams, q0, v0, tspan, dt; saveat)
Functional wrapper for the ODE solver.
"""
function solve_fput(p::SystemParams, q0, v0, tspan, dt; saveat=nothing)
    k, m = make_system(p)
    # inv_m avoids division per site per evaluation (18 evaluations per step);
    # F is the force buffer per bond, unique to this call ⇒ thread-safe.
    inv_m = 1.0 ./ m
    F     = Vector{Float64}(undef, p.boundary == :fixed ? p.N + 1 : p.N)
    p_ode = (k, inv_m, p.alpha, p.beta, p.boundary, F)
    
    # DifferentialEquations expects (v0, q0) for SecondOrderODE
    prob = SecondOrderODEProblem(fput_forces!, v0, q0, tspan, p_ode)
    
    # We use KahanLi8 (Symplectic) as requested by the original high-precision setup
    sol = solve(prob, KahanLi8(), dt=dt, saveat=saveat)
    
    # Return split trajectories for easier analysis
    # sol.u[i].x[1] is velocity, .x[2] is position
    V = stack(u.x[1] for u in sol.u)
    Q = stack(u.x[2] for u in sol.u)
    
    return Q, V, sol.t, k, m
end

end # module
