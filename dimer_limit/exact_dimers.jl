"""
    exact_dimers.jl

Exact reference at Δκ = 1 (PBC): κ₁ = 2 on the odd bonds, κ₂ = 0 on the even ones,
so the chain is N/2 independent dimers. The modal pipeline cannot describe this
point — the acoustic band collapses to N/2 degenerate zero modes, so `scaled_t`
has no clock and the modal basis is arbitrary inside the degenerate subspace —
so the reference is built in dimer coordinates instead.

With unit masses, each dimer (q₁, q₂) separates exactly into

    Q = (q₁+q₂)/2      Q̈ = 0                         (free translation)
    r =  q₂−q₁         r̈ = −2κ₁ r (1 + αr + βr²)     (anharmonic oscillator, ω = 2 for small r)

with conserved dimer energy  h = ṙ²/4 + κ₁(r²/2 + αr³/3 + βr⁴/4)  (reduced mass ½).

The script integrates the *full* 2N-dimensional chain with the production
integrator (KahanLi8) and checks it against the reduced dynamics solved
independently (Vern9, tolerance 1e-13). Agreement means Δκ = 1 is exactly
decoupled in the code, so every dimer energy is conserved: no thermalization,
T_th = ∞. That is the endpoint the Δκ < 1 runs approach.

Usage:
  julia --project=. dimer_limit/exact_dimers.jl
"""

using DifferentialEquations, Random, Printf
const ROOT = joinpath(@__DIR__, "..")
include(joinpath(ROOT, "src/fput_core.jl"));        using .FPUTCore
include(joinpath(ROOT, "src/fput_fast_runner.jl")); using .FPUTFastRunner

const N, α, β = 64, 0.1, 0.0
const T, DT, NSAVE = 2000.0, 0.05, 400
const κ₁ = 2.0                     # 1 + Δκ at Δκ = 1
const ND = N ÷ 2                   # number of dimers

bond_potential_dimer(r) = κ₁ * (r^2 / 2 + α * r^3 / 3 + β * r^4 / 4)
dimer_energy(r, ṙ) = ṙ^2 / 4 + bond_potential_dimer(r)

reduced_rhs!(du, u, p, t) = (du[1] = u[2]; du[2] = -2κ₁ * u[1] * (1 + u[1] * (α + β * u[1])))

function main()
    rng = Xoshiro(1)
    r₀, ṙ₀ = 2 .* (rand(rng, ND) .- 0.5), 0.5 .* randn(rng, ND)     # α|r| ≲ 0.1: weakly nonlinear
    Q₀, Q̇₀ = randn(rng, ND), 0.1 .* randn(rng, ND)

    # Dimer j = masses (2j−1, 2j); interleave back into chain order.
    q0, v0 = zeros(N), zeros(N)
    q0[1:2:end], q0[2:2:end] = Q₀ .- r₀ ./ 2, Q₀ .+ r₀ ./ 2
    v0[1:2:end], v0[2:2:end] = Q̇₀ .- ṙ₀ ./ 2, Q̇₀ .+ ṙ₀ ./ 2

    # Full chain, production integrator, at exactly Δκ = 1.
    sp = SystemParams(N, 1.0, 0.0, α, β, :periodic)
    Q, V, t, _, _ = solve_fput(sp, q0, v0, (0.0, T), DT; saveat = range(0, T; length = NSAVE + 1))

    r_full  = Q[2:2:end, :] .- Q[1:2:end, :]
    ṙ_full  = V[2:2:end, :] .- V[1:2:end, :]
    Qc_full = (Q[2:2:end, :] .+ Q[1:2:end, :]) ./ 2

    # Independent reference: each dimer's reduced ODE, high-order adaptive; COM is ballistic.
    r_ref = similar(r_full)
    for j in 1:ND
        sol = solve(ODEProblem(reduced_rhs!, [r₀[j], ṙ₀[j]], (0.0, T)), Vern9();
                    abstol = 1e-13, reltol = 1e-13, saveat = t)
        r_ref[j, :] = sol[1, :]
    end
    Qc_ref = Q₀ .+ Q̇₀ .* t'

    h = dimer_energy.(r_full, ṙ_full)
    dev_r, dev_c = maximum(abs, r_full .- r_ref), maximum(abs, Qc_full .- Qc_ref)
    drift = maximum(abs.(h .- h[:, 1:1]) ./ h[:, 1:1])

    @printf("N=%d (%d dimers), α=%g, β=%g, T=%g, dt=%g\n", N, ND, α, β, T, DT)
    @printf("  max |r_full − r_reduced|          = %.2e   (relative dimer coordinate)\n", dev_r)
    @printf("  max |Q_full − Q_ballistic|        = %.2e   (centre of mass)\n", dev_c)
    @printf("  max relative drift of h_j(t)      = %.2e   (energy exchange between dimers)\n", drift)
    @printf("  dimer energies h_j span           = [%.3g, %.3g]\n", minimum(h[:, 1]), maximum(h[:, 1]))

    ok = dev_r < 1e-6 && dev_c < 1e-9 && drift < 1e-6
    println(ok ? "\nPASS: Δκ = 1 is exactly N/2 decoupled dimers; every h_j is conserved." :
                 "\nFAIL: the full chain departs from the decoupled-dimer dynamics.")
    ok || exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
