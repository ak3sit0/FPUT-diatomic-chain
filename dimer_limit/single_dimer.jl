"""
    single_dimer.jl

Δκ = 1 (κ₁ = 2, κ₂ = 0), PBC: the chain is N/2 independent dimers. Excite the
relative coordinate of dimer 1 only (zero centre-of-mass motion, every other mass
at rest), integrate to 10⁵ cycles with the production integrator, check stability
and energy conservation, and compute the spectral entropy S̄(t).

Clock: there is no acoustic band at Δκ = 1, so "cycles" are periods of the dimer
oscillation, ω_d = √(2κ₁) = 2 (the flat optical band): t_cycles = t·ω_d/2π.

Spectral basis. Both bands are flat (ω = 0 and ω = 2, each N/2-fold degenerate),
so the eigenvectors `find_normal_modes` returns are an arbitrary rotation inside
each degenerate subspace — an entropy on them measures LAPACK's choice, not the
physics. The physical basis is the dimer basis (relative + centre of mass of each
dimer). S̄ is computed in both, and the difference is reported.

Usage:
  julia --project=. dimer_limit/single_dimer.jl [N=64] [E=0.445] [cycles=1e5]
"""

using LinearAlgebra, JLD2, Printf, LaTeXStrings
const ROOT = joinpath(@__DIR__, "..")
include(joinpath(ROOT, "src/fput_core.jl"));        using .FPUTCore
include(joinpath(ROOT, "src/fput_fast_runner.jl")); using .FPUTFastRunner
include(joinpath(ROOT, "src/fput_analysis.jl"));    using .FPUTAnalysis
include(joinpath(ROOT, "src/plotting_utils.jl"));   using .PlottingUtils
include(joinpath(ROOT, "src/plot_style.jl"));       using .PlotStyle

arg(i, d) = length(ARGS) >= i ? parse(Float64, ARGS[i]) : d

const α, β, DT = 0.1, 0.0, 0.05
const κ₁ = 2.0
const ω_d = sqrt(2κ₁)                 # small-amplitude dimer frequency
const NSAVE = 50_000                  # uniform: spectral_entropy's window is index-based

"""Total Hamiltonian of the periodic chain (unit masses)."""
function hamiltonian(q, v, k)
    N = length(q)
    0.5 * sum(abs2, v) + sum(bond_potential(k[i], q[mod1(i+1, N)] - q[i], α, β) for i in 1:N)
end

"""
Energies in the dimer basis, N × nt: row 2j−1 is the relative-coordinate energy of
dimer j (harmonic part, the analogue of a modal energy), row 2j its centre-of-mass
kinetic energy.
"""
function dimer_energies(Q, V)
    r, ṙ = Q[2:2:end, :] .- Q[1:2:end, :], V[2:2:end, :] .- V[1:2:end, :]
    Q̇    = (V[2:2:end, :] .+ V[1:2:end, :]) ./ 2
    E = similar(Q)
    E[1:2:end, :] .= ṙ .^ 2 ./ 4 .+ κ₁ .* r .^ 2 ./ 2      # reduced mass ½
    E[2:2:end, :] .= Q̇ .^ 2                               # total mass 2
    E
end

function main()
    N, E_in, cycles = Int(arg(1, 64)), arg(2, 0.445), arg(3, 1e5)
    T = cycles * 2π / ω_d

    # Dimer 1 = masses (1, 2), joined by the κ₁ bond. Equal and opposite displacements:
    # pure relative motion, zero centre of mass. Harmonic amplitude for energy E0.
    r0 = sqrt(2 * E_in / κ₁)          # NB: `2E0` would parse as the literal 2e0
    q0, v0 = zeros(N), zeros(N)
    q0[1], q0[2] = -r0 / 2, r0 / 2

    sp   = SystemParams(N, 1.0, 0.0, α, β, :periodic)
    k, m = make_system(sp)
    @printf("Δκ=1, PBC, N=%d, α=%g: dimer 1 with E=%g (r₀=%.3f, α·r₀=%.3f), T=%.3e (%.0e cycles of ω_d=2), dt=%g\n",
            N, α, E_in, r0, α*r0, T, cycles, DT)
    @printf("bonds: k[1]=%g (dimer 1), k[2]=%g (link to dimer 2)\n", k[1], k[2])

    t_wall = @elapsed Q, V, t, _, _ = solve_fput(sp, q0, v0, (0.0, T), DT;
                                                saveat = range(0, T; length = NSAVE + 1))
    @printf("integrated in %.1f s\n\n", t_wall)

    # ── Checks ────────────────────────────────────────────────────────────────
    H  = [hamiltonian(Q[:, i], V[:, i], k) for i in axes(Q, 2)]
    rel = abs.(H .- H[1]) ./ H[1]
    dH = maximum(rel)
    # Symplectic signature: the error oscillates at a bounded level instead of
    # growing. Compare the first and last tenth of the run.
    dec = length(rel) ÷ 10
    dH_first, dH_last = maximum(rel[2:dec]), maximum(rel[end-dec+1:end])
    finite  = all(isfinite, Q) && all(isfinite, V)
    others  = maximum(abs, Q[3:end, :]), maximum(abs, V[3:end, :])
    com     = maximum(abs, Q[1, :] .+ Q[2, :]), maximum(abs, V[1, :] .+ V[2, :])
    r       = Q[2, :] .- Q[1, :]

    @printf("finite everywhere                : %s\n", finite)
    @printf("energy H(0)                      : %.12f\n", H[1])
    @printf("max relative energy error |ΔH|/H : %.2e\n", dH)
    @printf("  first 10%% vs last 10%% of run   : %.2e vs %.2e  (bounded ⇔ comparable)\n", dH_first, dH_last)
    @printf("max |q|,|v| on masses 3…N        : %.2e, %.2e  (must stay 0: κ₂ = 0)\n", others...)
    @printf("max |q₁+q₂|, |v₁+v₂| (dimer COM) : %.2e, %.2e\n", com...)
    @printf("relative coordinate r ∈          : [%.4f, %.4f]  (α·r turnover at %.0f)\n",
            minimum(r), maximum(r), -1/α)

    ok = finite && dH < 1e-5 && dH_last < 3 * dH_first &&
         maximum(others) < 1e-12 && maximum(com) < 1e-9
    println(ok ? "\nStable, energy conserved, energy confined to dimer 1.\n" :
                 "\nCHECK FAILED — not computing the entropy.\n")
    ok || exit(1)

    # ── Spectral entropy ──────────────────────────────────────────────────────
    t_cyc = t .* ω_d ./ 2π
    S_dimer = spectral_entropy(dimer_energies(Q, V), SMOOTH_DELTA)

    freq, Vm = find_normal_modes(k, m, :periodic)
    sm       = sqrt.(m)   # compute_modal_energies on all columns at once
    E_modal  = 0.5 .* ((Vm' * (sm .* V)) .^ 2 .+ freq .^ 2 .* (Vm' * (sm .* Q)) .^ 2)
    S_modal  = spectral_entropy(E_modal, SMOOTH_DELTA)

    @printf("S̄ dimer basis   : initial %.3e, final %.3e, max %.3e   (ln N = %.3f)\n",
            S_dimer[2], S_dimer[end], maximum(S_dimer), log(N))
    @printf("S̄ LAPACK modes  : initial %.3f, final %.3f, max %.3f   (bands are degenerate, so this basis is not guaranteed; here it agrees)\n",
            S_modal[2], S_modal[end], maximum(S_modal))
    @printf("degeneracy      : %d modes at ω≈0, %d at ω≈2\n",
            count(<(1e-6), freq), count(f -> abs(f - 2) < 1e-6, freq))

    outdir = joinpath(@__DIR__, "data"); mkpath(outdir)
    jldsave(joinpath(outdir, "single_dimer_N$(N).jld2");
            t_cycles = t_cyc, H = H, r = r, S_dimer = S_dimer, S_modal = S_modal,
            N = N, E_in = E_in, alpha = α, dt = DT)

    # ── Figure ────────────────────────────────────────────────────────────────
    apply_style!()
    keep = t_cyc .> 0
    curves = [Curve(logdownsample(t_cyc[keep], S_dimer[keep])..., "dimer basis";
                    color = PALETTE_DELTA[1], linewidth = 2.5),
              Curve(logdownsample(t_cyc[keep], S_modal[keep])..., "normal-mode basis (degenerate)";
                    color = PALETTE_DELTA[5], linestyle = :dash, linewidth = 2.0)]
    fig = logplot(; xlabel = L"t\ \mathrm{(cycles\ of\ }\omega_d=2)", ylabel = L"\bar{S}(t)",
                  decades = 0:Int(log10(cycles)), legend = :right,
                  # Pinned to [0, ln N]: S̄ is ~1e-15 here, and autoscaling would blow
                  # roundoff up into what looks like structure.
                  ylims = (-0.05 * log(N), 1.08 * log(N)),
                  title = latexstring("\\Delta\\kappa = 1,\\ N = $N,\\ \\mathrm{one\\ dimer\\ excited}"))
    draw!(fig, curves)
    guide_hline!(fig, log(N); annotation = L"\ln N\ \mathit{(equipartition)}", side = :left)
    save_fig(fig, joinpath(@__DIR__, "figures"), "single_dimer_entropy_N$(N)"; exts = (".png",))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
