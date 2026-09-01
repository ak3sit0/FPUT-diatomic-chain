"""
    compute_ftmle.jl

Finite-time maximum Lyapunov exponent via the two-trajectory method of Benettin
et al.: integrate a reference and a perturbed trajectory in chunks of `T_renorm`,
renormalizing the separation back to `delta0` each time.

Unlike the sweep scripts this writes **one JLD2 per Δκ** (with `t_cycles`,
`lambda`, `omega_ref`, …), not one collective `results` file, so it does not use
`SweepDriver`. It does share the config and case setup, which is where its own
copies of `read_config` and `mode_ic` used to live.

Usage:
  julia -t <n> --project=. scripts/compute/compute_ftmle.jl [config.toml] [T_max] [T_renorm] [d1,d2,...]
"""

using JLD2, LinearAlgebra, Random, Dates, Base.Threads
include("../../src/fput_core.jl");        using .FPUTCore
include("../../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../../src/experiment.jl");       using .Experiment
include("../../src/case_setup.jl");       using .CaseSetup

const DEFAULT_CONFIG = "configs/production/periodic_N64.toml"

"""
    run_ftmle(case, q_ref, v_ref, dt; T_max, T_renorm, delta0, rng_seed) -> (t, λ)

Benettin protocol. λ(t) is the running average `Λ/t` of the logarithmic
separation growth, so it decays as `t^{-1}` for regular motion and plateaus at
the Lyapunov exponent for chaotic motion.
"""
function run_ftmle(case::Case, q_ref, v_ref, dt::Float64;
                   T_max = 1e7, T_renorm = 200.0, delta0 = 1e-8, rng_seed = 1234)
    N  = case.N
    sp = case.sp

    # Random unit perturbation in the 2N-dimensional phase space
    rng = MersenneTwister(rng_seed)
    dv  = randn(rng, 2N)
    dv ./= norm(dv)
    q_per = q_ref .+ delta0 .* dv[1:N]
    v_per = v_ref .+ delta0 .* dv[N+1:end]

    n_chunks = Int(ceil(T_max / T_renorm))
    t_vec    = Vector{Float64}(undef, n_chunks)
    λ_vec    = Vector{Float64}(undef, n_chunks)
    Λ, t     = 0.0, 0.0

    println("ftMLE: N=$N  Δκ=$(case.delta)  mode=$(case.ref_mode)  T_max=$T_max  chunks=$n_chunks")
    flush(stdout)

    for i in 1:n_chunks
        t_next = t + T_renorm
        Q1, V1, _, _, _ = solve_fput(sp, q_ref, v_ref, (t, t_next), dt; saveat = [t_next])
        Q2, V2, _, _, _ = solve_fput(sp, q_per, v_per, (t, t_next), dt; saveat = [t_next])

        q_ref, v_ref = Q1[:, end], V1[:, end]
        Δq = Q2[:, end] .- q_ref
        Δv = V2[:, end] .- v_ref
        d  = norm(vcat(Δq, Δv))

        Λ += log(d / delta0)
        t  = t_next
        t_vec[i], λ_vec[i] = t, Λ / t

        # Renormalize: bring the perturbed trajectory back to distance delta0
        q_per = q_ref .+ (delta0 / d) .* Δq
        v_per = v_ref .+ (delta0 / d) .* Δv

        if mod(i, max(1, n_chunks ÷ 10)) == 0
            println("  chunk $i/$n_chunks  t=$(round(t; digits=1))  λ=$(round(Λ/t; sigdigits=4))")
            flush(stdout)
        end
    end

    t_vec, λ_vec
end

function save_result(spec, case, config_path, T_max, t_vec, λ_vec)
    outdir = joinpath(spec.base_dir, "ftmle")
    mkpath(outdir)
    tag  = replace(string(case.delta), "." => "p")
    path = joinpath(outdir, "ftmle_$(spec.boundary)_delta$(tag)_$(Dates.today()).jld2")
    jldsave(path;
        t_physical  = t_vec,
        t_cycles    = t_vec .* case.omega_ref ./ (2π),
        lambda      = λ_vec,
        omega_ref   = case.omega_ref,
        delta_k     = case.delta,
        boundary    = string(spec.boundary),
        init_mode   = case.ref_mode,
        T_max       = T_max,
        config_path = config_path,
    )
    println("Saved: $path")
end

function main()
    config   = get(ARGS, 1, DEFAULT_CONFIG)
    T_max    = parse(Float64, get(ARGS, 2, "10000000.0"))
    T_renorm = parse(Float64, get(ARGS, 3, "200.0"))

    spec = parse_spec(config)

    deltas = if length(ARGS) >= 4
        requested = parse.(Float64, split(ARGS[4], ","))
        for d in requested
            any(x -> isapprox(x, d; atol = 1e-12), spec.delta_values) ||
                error("delta=$d is not in delta_values=$(spec.delta_values)")
        end
        requested
    else
        spec.delta_values
    end

    N     = first(spec.N_values)
    param = first(spec.param_values)
    println("Cases: Δκ = $deltas   N=$N   ($(nthreads()) threads)")
    flush(stdout)

    @threads for Δ in deltas
        case   = build_case(spec, N, param, Δ)
        q0, v0 = initial_condition(case, spec, spec.seed_base)
        t_vec, λ_vec = run_ftmle(case, q0, v0, spec.DT;
                                 T_max = T_max, T_renorm = T_renorm)
        save_result(spec, case, config, T_max, t_vec, λ_vec)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
