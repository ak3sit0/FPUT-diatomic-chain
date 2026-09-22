# Energy-error check for the robustness runs: N=256, PBC, DT=0.2, KahanLi8.
# Symplectic ⇒ bounded oscillating error, so ~1e4 cycles characterizes the bound.
# Usage: julia --project=. scripts/compute/check_energy_drift.jl [N] [Δκ] [cycles]
include("../../src/fput_core.jl");        using .FPUTCore
include("../../src/fput_fast_runner.jl"); using .FPUTFastRunner
using LinearAlgebra, Printf

N      = length(ARGS) ≥ 1 ? parse(Int, ARGS[1])       : 256
Δκ     = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2])   : 0.7
cycles = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3])   : 1.0e4
α, ε, DT, mode = 0.1, 0.00695, 0.2, 2

sp   = SystemParams(N, Δκ, 0.0, α, 0.0, :periodic)
k, m = make_system(sp)
freq, V = find_normal_modes(k, m, :periodic)
E0   = ε * N
q0   = sqrt(2E0) / freq[mode] .* (Diagonal(1 ./ sqrt.(m)) * V)[:, mode]
v0   = zeros(N)

energy(q, v) = sum(@. m * v^2 / 2) +
               sum(bond_potential(k[i], q[mod1(i + 1, N)] - q[i], α, 0.0) for i in 1:N)

TMAX = cycles * 2π / freq[mode]
H0   = energy(q0, v0)
Q, Vv, t, _, _ = solve_fput(sp, q0, v0, (0.0, TMAX), DT; saveat = range(0, TMAX; length = 2001))
relerr = [abs(energy(Q[:, j], Vv[:, j]) - H0) / H0 for j in eachindex(t)]

@printf("N=%d Δκ=%.2f DT=%.2f cycles=%.1e  ω_max·DT=%.3f\n", N, Δκ, DT, cycles, maximum(freq) * DT)
@printf("max|ΔE/E| = %.3e   |ΔE/E| at end = %.3e\n", maximum(relerr), relerr[end])
@printf("max over first half = %.3e, second half = %.3e  (no drift ⇔ similar)\n",
        maximum(relerr[1:end÷2]), maximum(relerr[end÷2:end]))
