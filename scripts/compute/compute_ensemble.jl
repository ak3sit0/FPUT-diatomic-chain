"""
    compute_ensemble.jl

Phase-ensemble sweep over `N × param × Δ`: `n_real` realizations per case with
random phases on a selected band, averaged, plus thermalization statistics.

Replaces `compute_ensemble.jl` + `compute_ensemble_Nsweep.jl`. The two differed
in where the parallelism lives, and that is now derived from `n_real` rather than
being a separate script (see `main`). The N-sweep features — `scaled_t_max` and
`n_samples` resolved per N — are handled generically by `Experiment.resolve_budget`,
and `T_therm_*` is now emitted in both cases, closing the schema mismatch noted in
`docs/physics_diagnostics.md`.

Energy convention:
  `energy_density` ⇒ E_total = ε·N  (controls the sweep; preferred)
  `initial_energy` ⇒ E_total fixed, so ε falls as 1/N  (legacy)

Usage:
  julia --project=. scripts/compute/compute_ensemble.jl <config.toml>
"""

using JLD2, Statistics, LinearAlgebra, Base.Threads
include("../../src/fput_core.jl");        using .FPUTCore
include("../../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../../src/fput_analysis.jl");    using .FPUTAnalysis
include("../../src/block_integration.jl");using .BlockIntegration
include("../../src/experiment.jl");       using .Experiment
include("../../src/case_setup.jl");       using .CaseSetup
include("../../src/sweep_driver.jl");     using .SweepDriver

const USAGE = "Usage: julia --project=. scripts/compute/compute_ensemble.jl <config.toml>"

# ── One realization ───────────────────────────────────────────────────────────

"""
    run_realization(case, spec, seed, label) -> NamedTuple | nothing

One trajectory, tracking per-band energies through `on_block`. Returns `nothing`
when the run goes unstable (total energy >100× its initial value, or non-finite),
so the caller can drop it instead of poisoning the ensemble average.
"""
function run_realization(case::Case, spec, seed::Integer, label::AbstractString)
    E_ac_blocks, E_opt_blocks = Vector{Vector{Float64}}(), Vector{Vector{Float64}}()

    # idx_ds selects the (downsampled) saved columns within this block's modal_Eb;
    # summing over case.k_ac/k_opt rows collapses per-mode energy to per-band energy.
    on_block = (modal_Eb, idx_ds) -> begin
        push!(E_ac_blocks, vec(sum(modal_Eb[case.k_ac, idx_ds], dims = 1)))
        isempty(case.k_opt) || push!(E_opt_blocks, vec(sum(modal_Eb[case.k_opt, idx_ds], dims = 1)))
    end

    check_instability = modal_Eb -> begin
        E_cur = sum(@view modal_Eb[:, end])  # total energy at the block's last saved column
        bad = E_cur > 100 * case.budget.E_total || !isfinite(E_cur)
        bad && println("  $label INSTABILITY (E=$(round(E_cur; sigdigits=3)) vs E_total=$(case.budget.E_total)); aborting")
        bad
    end

    q0, v0 = initial_condition(case, spec, seed)
    b = case.budget
    r = integrate_in_blocks(case.sp, q0, v0, case.N, case.freq, case.V, case.m,
                            case.ref_mode;
                            TMAX = b.TMAX, T_block = b.T_block, DT = spec.DT,
                            save_every = b.save_every, downsample = spec.downsample,
                            track_abs_time = true, debug = spec.debug, label = label,
                            on_block = on_block, check_instability = check_instability)

    r.aborted && return nothing

    (; scaled_t = r.scaled_t, t_abs = r.t_abs, modal_E = r.modal_E,
       entropy = FPUTAnalysis.spectral_entropy(r.modal_E, spec.entropy_delta),
       E_acoustic = isempty(E_ac_blocks)  ? Float64[] : reduce(vcat, E_ac_blocks),
       E_optical  = isempty(E_opt_blocks) ? Float64[] : reduce(vcat, E_opt_blocks))
end

# ── Thermalization time ───────────────────────────────────────────────────────

"""
    thermalization_stats(entropies, scaled_t, k_band, N; xi_th = 0.5) -> NamedTuple

Per-realization crossing time of the effective mode occupancy ξ(t) = e^{S(t)}/N,
then ensemble statistics. `T_th` is the first time ξ exceeds `xi_th`.

Huang et al. (2016) / Wang et al. (2020, 2024) convention: equipartition is
ξ_eq = 1/2, so the threshold is the *absolute* entropy `S > log(N·xi_th) = log(N/2)`,
not a fractional rise from the initial band entropy. This replaces the earlier
`f = 1 - 1/e` rise on `(S - S₀)/(log N - S₀)`, which mixed the band size into the
threshold — ξ₀ depends on |k_band| but ξ_eq does not, so the absolute form is what
makes runs with different band sizes comparable.

Note this is *not* the same quantity as `PlottingUtils.thermalization_time`,
which thresholds the optical energy instead — see `docs/physics_diagnostics.md`.
"""
function thermalization_stats(entropies, scaled_t, k_band, N; xi_th::Real = 0.5)
    S_th = log(N * xi_th)  # ξ > xi_th ⟺ S > log(N·xi_th); for xi_th=1/2 this is log(N/2)

    T_vec = map(entropies) do S_r
        nt  = min(length(S_r), length(scaled_t))  # guard against off-by-one length mismatches
        idx = findfirst(>(S_th), @view S_r[1:nt])  # first sample with ξ above threshold
        isnothing(idx) ? Inf : scaled_t[idx]       # never thermalized within the run ⇒ Inf
    end

    finite = filter(isfinite, T_vec)
    (; T_therm_mean   = isempty(finite) ? Inf : mean(finite),
       T_therm_std    = length(finite) > 1 ? std(finite) : NaN,
       T_therm_median = isempty(finite) ? Inf : median(finite),
       T_therm_vec    = T_vec,
       frac_therm     = length(finite) / length(T_vec),
       n_therm        = length(finite),
       threshold_xi   = xi_th,
       xi_initial     = length(k_band) / N)  # ξ₀, for reporting: varies with band size
end

# ── One case (the whole ensemble) ─────────────────────────────────────────────

"""
    run_ensemble_case(spec, task, inner_parallel) -> NamedTuple | nothing

Run `n_real` realizations and reduce them.

Realizations are launched in batches of `nthreads()` rather than all at once:
with `n_real=100` and a 13 MB `modal_E`, launching all of them would keep ~1.3 GB
alive per case. The reduction is serial and in order of `r`, so the result does
not depend on which thread finishes first — that is what makes a seeded run
reproducible bit for bit.
"""
function run_ensemble_case(spec, task, inner_parallel::Bool)
    label_case = "N=$(task.N) p=$(task.param) Δ=$(task.delta)"
    case = build_case(spec, task.N, task.param, task.delta)
    b    = case.budget
    println("[$label_case] ε=$(round(b.E_total / case.N; sigdigits=4))  " *
            "TMAX=$(round(b.TMAX; sigdigits=4))  save_every=$(b.save_every)  " *
            "n_real=$(spec.n_real)  band=$(length(case.excited)) modes")

    entropies   = Vector{Vector{Float64}}()
    E_opt_reals = Vector{Vector{Float64}}()
    modal_acc = E_ac_acc = E_opt_acc = nothing
    t_ref = nothing
    n_ok  = 0

    batch = inner_parallel ? max(nthreads(), 1) : 1  # batch=1 ⇒ effectively a plain serial loop
    for start in 1:batch:spec.n_real
        rs    = start:min(start + batch - 1, spec.n_real)  # last batch may be shorter than nthreads()
        tasks = Vector{Any}(undef, length(rs))

        for (i, r) in enumerate(rs)
            lbl = "[$label_case r=$r/$(spec.n_real)]"
            # spec.seed_base + r (not i or start+i) keeps each realization's seed tied to its
            # global index r, so results are reproducible regardless of batch size/thread count.
            tasks[i] = inner_parallel ?
                Threads.@spawn(run_realization(case, spec, spec.seed_base + r, lbl)) :
                run_realization(case, spec, spec.seed_base + r, lbl)
        end

        for (i, r) in enumerate(rs)
            res = inner_parallel ? fetch(tasks[i]) : tasks[i]
            tasks[i] = nothing              # free as soon as possible

            isnothing(res) && (println("  [r=$r] discarded (unstable)"); continue)

            n_ok += 1
            push!(entropies, res.entropy)
            push!(E_opt_reals, res.E_optical)
            t_ref = res.scaled_t

            if isnothing(modal_acc)
                # copy(): res.modal_E is owned by this realization's integration buffer and
                # would otherwise get mutated/reused by a later thread.
                modal_acc, E_ac_acc, E_opt_acc =
                    copy(res.modal_E), copy(res.E_acoustic), copy(res.E_optical)
            else
                # An unstable/aborted realization can save fewer blocks than its predecessors,
                # so accumulate only over the common (shortest) length seen so far.
                nt = min(size(res.modal_E, 2), size(modal_acc, 2))
                modal_acc = modal_acc[:, 1:nt] .+ res.modal_E[:, 1:nt]
                E_ac_acc  = _accumulate(E_ac_acc, res.E_acoustic)
                E_opt_acc = _accumulate(E_opt_acc, res.E_optical)
            end
        end
    end

    n_ok == 0 && (println("[$label_case] all realizations unstable; case discarded"); return nothing)
    n_ok < spec.n_real && println("  Warning: only $n_ok/$(spec.n_real) realizations stable")

    modal_acc ./= n_ok  # elementwise mean over the n_ok stable realizations (division, not sum)
    E_ac_acc  ./= n_ok
    E_opt_acc  ./= n_ok

    nt_min  = minimum(length, entropies)                        # shortest surviving trajectory
    ent_mat = reduce(hcat, [e[1:nt_min] for e in entropies])'    # hcat stacks as columns, so
    S_mean  = vec(mean(ent_mat, dims = 1))                       # the trailing ' makes rows =
    S_std   = vec(std(ent_mat, dims = 1))                        # realizations, dims=1 reduces those

    th = thermalization_stats(entropies, t_ref, case.excited, case.N)

    # ξ̄ rather than S/logN: ξ = e^S/N is the quantity the ξ>1/2 threshold is stated in,
    # so the printed value can be read directly against it.
    println("[$label_case] S̄=$(round(S_mean[end]; digits=4)) ± $(round(S_std[end]; digits=4))  " *
            "ξ̄=$(round(exp(S_mean[end])/case.N; digits=4))  " *
            "T_th(med)=$(round(th.T_therm_median; sigdigits=3))  therm=$(th.n_therm)/$(spec.n_real)")

    (; N = case.N, param = case.param, Delta = case.delta,
       scaled_t = t_ref,
       entropy_mean = S_mean, entropy_std = S_std,
       modal_E_mean = modal_acc,
       E_acoustic_mean = E_ac_acc, E_optical_mean = E_opt_acc,
       entropy_realizations = entropies, E_optical_realizations = E_opt_reals,
       th...,  # splices T_therm_mean/std/median/vec, frac_therm, n_therm, threshold_xi, xi_initial
       n_real = spec.n_real, seed_base = spec.seed_base,
       branch = spec.branch, k_band = case.excited,
       E_total = b.E_total, energy_density = b.E_total / case.N,
       omega_ref = case.omega_ref, TMAX = b.TMAX, save_every = b.save_every,
       init_type = spec.init_type)
end

# Sums, truncating to the shorter of the two (see the aborted-realization note above);
# either side can be empty when a case has no optical band (fixed boundary), so pass acc through.
_accumulate(acc, new) = isempty(acc) || isempty(new) ? acc :
                        (nt = min(length(acc), length(new)); acc[1:nt] .+ new[1:nt])

# ── Main ──────────────────────────────────────────────────────────────────────

function main()
    isempty(ARGS) && error(USAGE)
    spec  = parse_spec(ARGS[1])
    tasks = sweep_tasks(spec)

    # Where the parallelism goes depends on n_real, which is what used to separate
    # this script from its N-sweep twin:
    #   n_real == 1  → nothing to parallelize inside, so parallelize over cases.
    #   n_real  > 1  → the inner loop already saturates the node; nesting @threads
    #                  outside would just block threads waiting.
    inner = spec.n_real > 1
    outer = inner ? :serial : :dynamic

    println("$(length(tasks)) cases, n_real=$(spec.n_real), $(nthreads()) threads — " *
            "parallelism " * (inner ? "over realizations (cases in series)" :
                                      "over cases (one realization each)") * "\n")

    run_sweep(t -> run_ensemble_case(spec, t, inner), tasks;
              outfile     = output_path(spec, "ensemble_results"),
              config_path = ARGS[1],
              parallel    = outer)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
