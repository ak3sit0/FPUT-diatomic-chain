"""
    spectrum.jl

Preflight for the Δκ → 1 study. For each Δκ, prints what the pipeline will face
before any CPU is spent: band edges, the acoustic clock ω₂, the time budget, and
how close the initial condition sits to the turnover of the cubic potential.

Everything is read off `find_normal_modes`, not from closed forms, so it stays
correct if the chain definition changes. (The closed forms — ω_ac,max = √(2(1−Δκ)),
ω_op,min = √(2(1+Δκ)), N·ω₂ → 2π√(1−Δκ²) — were checked symbolically and are
reproduced by the table.)

Usage:
  julia --project=. dimer_limit/spectrum.jl [N=64] [alpha=0.1] [E_total=0.445] [scaled_t_max=1e3] [target_strain=0.2]
"""

using LinearAlgebra, Printf
const ROOT = joinpath(@__DIR__, "..")
include(joinpath(ROOT, "src/fput_core.jl")); using .FPUTCore

const DELTAS = [0.5, 0.9, 0.99, 0.999, 0.9999, 1.0]
const DT = 0.05
const SEC_PER_STEP_SITE = 105e-9   # FPUT cost model: wall ≈ 105 ns · N · TMAX/dt

arg(i, default) = length(ARGS) >= i ? parse(Float64, ARGS[i]) : default

"""
    max_strain(N, α, E, freq, V, m) -> Float64

Largest `α·|d|` over all bonds at t = 0 for mode 2 carrying `E` (harmonic
amplitude √(2E)/ω₂). The cubic bond potential turns over at `α·d = −1`, so this
is the distance to escape; it grows as ∝ 1/√κ₂ at fixed E, because the acoustic
mode moves the dimers rigidly and stores all its energy in the soft bonds.
"""
function max_strain(N, α, E, freq, V, m)
    U = Diagonal(1 ./ sqrt.(m)) * V
    q = sqrt(2E) / freq[2] .* U[:, 2]
    α * maximum(abs(q[mod1(i+1, N)] - q[i]) for i in 1:N)
end

function human(seconds)
    seconds < 90    && return @sprintf("%.0f s", seconds)
    seconds < 5400  && return @sprintf("%.0f min", seconds / 60)
    seconds < 172800 && return @sprintf("%.1f h", seconds / 3600)
    @sprintf("%.0f d", seconds / 86400)
end

flag(s) = s < 0.3 ? "ok" : s < 0.6 ? "marginal" : "UNSAFE"

function main()
    N, α, E = Int(arg(1, 64)), arg(2, 0.1), arg(3, 0.445)
    st, target = arg(4, 1e3), arg(5, 0.2)

    @printf("PBC, N=%d, α=%g, mode 2 with E_total=%g (ε=%.3g), scaled_t_max=%g\n\n",
            N, α, E, E/N, st)
    @printf("%-8s %-9s %-9s %-8s %-10s %-9s %-9s %-9s %s\n",
            "Δκ", "ω_ac,max", "ω_op,min", "N·ω₂", "TMAX", "wall", "α·|d|max", "verdict",
            "E for α·|d|=$target")

    for dk in DELTAS
        sp = SystemParams(N, dk, 0.0, α, 0.0, :periodic)
        k, m = make_system(sp)
        freq, V = find_normal_modes(k, m, :periodic)

        if freq[2] < OMEGA_MIN
            @printf("%-8g %-9.3g %-9.3g %-8s %s\n", dk, freq[N÷2], freq[N÷2+1], "—",
                    "degenerate: $(count(<(OMEGA_MIN), freq)) modes at ω≈0, no acoustic clock " *
                    "(ref_frequency refuses) → dimer_limit/exact_dimers.jl")
            continue
        end

        TMAX = st * 2π / freq[2]
        s    = max_strain(N, α, E, freq, V, m)
        @printf("%-8g %-9.4g %-9.4g %-8.4g %-10.3g %-9s %-9.3g %-9s %.3g\n",
                dk, freq[N÷2], freq[N÷2+1], N*freq[2], TMAX,
                human(SEC_PER_STEP_SITE * N * TMAX / DT), s, flag(s), E * (target / s)^2)
    end

    println("\nverdict: ok < 0.3 ≤ marginal < 0.6 ≤ UNSAFE (heuristic thresholds on α·|d|max; " *
            "1.0 is the turnover of the cubic potential).")
    println("The last column is E_total that puts mode 2 at the target strain — the fixed-strain " *
            "alternative to fixed ε.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# Suggested test (add to test/runtests.jl):
# @testset "band edges match closed forms" begin
#     N = 64
#     for dk in (0.5, 0.9, 0.99)
#         freq, _ = find_normal_modes(make_system(SystemParams(N, dk, 0.0, 0.1, 0.0, :periodic))..., :periodic)
#         @test freq[N÷2]   ≈ sqrt(2(1 - dk)) rtol = 1e-8
#         @test freq[N÷2+1] ≈ sqrt(2(1 + dk)) rtol = 1e-8
#         @test N * freq[2] ≈ 2π * sqrt(1 - dk^2) rtol = 5e-3   # → equality as N → ∞
#     end
# end
