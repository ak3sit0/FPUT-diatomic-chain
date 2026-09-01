"""
    test_compute_pipeline.jl

Unit tests for the config → case → sweep pipeline that the compute scripts share
(`src/experiment.jl`, `src/case_setup.jl`, `src/sweep_driver.jl`), plus the two
guards added to `src/fput_core.jl`.

Written against the failure modes that actually bit this project:

- `init_mode` being read by one script and ignored by another (FBC ran mode 1
  while ftMLE ran mode 2 from the same TOML).
- A near-zero reference frequency going unnoticed because it is ~3e-8 rather
  than exactly 0, which collapses `scaled_t` and inflates `TMAX` by ~10⁶.
- Per-N quantities resolved once at parse time, which silently breaks any sweep
  where N varies.
"""

using Test, TOML, JLD2, LinearAlgebra

@testset "Compute pipeline" begin

# ── Helper: write a throwaway TOML ────────────────────────────────────────────

function write_toml(dir, name; physics = Dict(), simulation = Dict(), output = Dict())
    phys = merge(Dict("N" => 32, "boundary" => "periodic", "system_type" => "springs",
                      "nonlinear" => "alpha", "param_values" => [0.1],
                      "delta_values" => [0.1]), physics)
    sim  = merge(Dict("TMAX" => 100.0, "T_block" => 50.0, "DT" => 0.05,
                      "save_every" => 10, "downsample" => 1, "debug" => false), simulation)
    out  = merge(Dict("base_dir" => joinpath(dir, "data")), output)
    path = joinpath(dir, name)
    open(path, "w") do io
        TOML.print(io, Dict("physics" => phys, "simulation" => sim, "output" => out))
    end
    path
end

# ── FPUTCore guards ───────────────────────────────────────────────────────────

@testset "ref_frequency rejects the PBC zero mode" begin
    # Mode 1 under PBC is the uniform k=0 translation. Its eigenvalue is roundoff,
    # so ω comes out ~3e-8: finite, so nothing throws on its own, which is exactly
    # why this needs an explicit guard.
    N  = 64
    sp = SystemParams(N, 0.1, 0.0, 0.1, 0.0, :periodic)
    k, m    = make_system(sp)
    freq, _ = find_normal_modes(k, m, :periodic)

    @test freq[1] < 1e-6            # numerically zero…
    @test freq[1] != 0.0            # …but not exactly zero: the silent-failure trap
    @test_throws ErrorException ref_frequency(freq, 1)

    # Mode 2 is the lowest physical mode and must pass through untouched.
    @test ref_frequency(freq, 2) == freq[2]
    @test ref_frequency(freq, 2) > 1e-6

    # Out-of-range indices are caught rather than reaching freq[] as a BoundsError.
    @test_throws ErrorException ref_frequency(freq, 0)
    @test_throws ErrorException ref_frequency(freq, N + 1)
end

@testset "ref_frequency accepts mode 1 under fixed BC" begin
    # Fixed ends have no zero mode, so mode 1 is a legitimate lowest acoustic mode.
    sp = SystemParams(64, 0.1, 0.0, 0.1, 0.0, :fixed)
    k, m    = make_system(sp)
    freq, _ = find_normal_modes(k, m, :fixed)

    @test freq[1] > 1e-6
    @test ref_frequency(freq, 1) == freq[1]
end

@testset "find_normal_modes returns ascending frequencies" begin
    # Five scripts used to run sortperm(freq) after this call. LAPACK's symmetric
    # solver already returns ascending eigenvalues and sqrt is monotone, so those
    # sorts were no-ops — this pins the property the removal relies on.
    for bc in (:fixed, :periodic), dk in (0.0, 0.1, 0.5)
        sp = SystemParams(32, dk, 0.0, 0.1, 0.0, bc)
        k, m    = make_system(sp)
        freq, _ = find_normal_modes(k, m, bc)
        @test issorted(freq)
    end
end

@testset "derive_band_indices validates manual band limits" begin
    N  = 32
    sp = SystemParams(N, 0.1, 0.0, 0.1, 0.0, :periodic)
    k, m    = make_system(sp)
    freq, _ = find_normal_modes(k, m, :periodic)

    # The automatic branch deliberately starts at 2 under PBC; the manual branch
    # used to bypass that protection entirely.
    @test_throws ErrorException derive_band_indices("acoustic", N, freq, :periodic;
                                                    k_band_start = 1, k_band_end = 8)
    @test_throws ErrorException derive_band_indices("acoustic", N, freq, :periodic;
                                                    k_band_start = 0, k_band_end = 8)
    @test_throws ErrorException derive_band_indices("acoustic", N, freq, :periodic;
                                                    k_band_start = 2, k_band_end = N + 5)
    @test_throws ErrorException derive_band_indices("acoustic", N, freq, :periodic;
                                                    k_band_start = 8, k_band_end = 2)

    # A valid manual band still works.
    @test derive_band_indices("acoustic", N, freq, :periodic;
                              k_band_start = 2, k_band_end = 8) == collect(2:8)

    # The automatic acoustic band never includes the translation mode.
    @test first(derive_band_indices("acoustic", N, freq, :periodic)) == 2
end

# ── Experiment: parsing and policy ────────────────────────────────────────────

@testset "excited mode default is physical, override is honoured" begin
    mktempdir() do dir
        # No init_mode → default per boundary.
        periodic = parse_spec(write_toml(dir, "p.toml"; physics = Dict("boundary" => "periodic")))
        fixed    = parse_spec(write_toml(dir, "f.toml"; physics = Dict("boundary" => "fixed")))
        @test default_excited_mode(:periodic) == 2
        @test default_excited_mode(:fixed)    == 1
        @test excited_mode(periodic) == 2
        @test excited_mode(fixed)    == 1

        # Explicit init_mode wins — this is the field compute_trajectories.jl used
        # to read and then ignore, which is how FBC data diverged from ftMLE data.
        override = parse_spec(write_toml(dir, "o.toml";
            physics = Dict("boundary" => "fixed", "init_mode" => 3)))
        @test excited_mode(override) == 3
    end
end

@testset "legacy initial_condition key is accepted with a warning" begin
    mktempdir() do dir
        path = write_toml(dir, "legacy.toml";
                          physics = Dict("initial_condition" => "mode"))
        spec = @test_logs (:warn,) match_mode = :any parse_spec(path)
        @test spec.init_type == "mode"
    end
end

@testset "N_values falls back to N" begin
    mktempdir() do dir
        @test parse_spec(write_toml(dir, "a.toml"; physics = Dict("N" => 48))).N_values == [48]
        @test parse_spec(write_toml(dir, "b.toml";
                physics = Dict("N_values" => [16, 32, 64]))).N_values == [16, 32, 64]
    end
end

# ── Experiment: per-N budget resolution ───────────────────────────────────────

@testset "energy_density fixes ε, initial_energy fixes E_total" begin
    mktempdir() do dir
        # ε fixed: E_total must scale with N, keeping the nonlinearity comparable.
        eps_spec = parse_spec(write_toml(dir, "eps.toml";
                                         physics = Dict("energy_density" => 0.00695)))
        @test resolve_budget(eps_spec, 64,  0.1).E_total ≈ 0.00695 * 64
        @test resolve_budget(eps_spec, 256, 0.1).E_total ≈ 0.00695 * 256

        # E_total fixed (legacy): ε then falls as 1/N.
        tot_spec = parse_spec(write_toml(dir, "tot.toml";
                                         physics = Dict("initial_energy" => 0.445)))
        @test resolve_budget(tot_spec, 64,  0.1).E_total ≈ 0.445
        @test resolve_budget(tot_spec, 256, 0.1).E_total ≈ 0.445
    end
end

@testset "scaled_t_max makes TMAX scale like 1/ω_ref" begin
    mktempdir() do dir
        spec = parse_spec(write_toml(dir, "st.toml";
                                     simulation = Dict("scaled_t_max" => 1000.0)))
        # ω_ref ≈ 6.25/N, so halving ω_ref (doubling N) doubles TMAX.
        b1 = resolve_budget(spec, 64,  0.0976)
        b2 = resolve_budget(spec, 128, 0.0488)
        @test b1.TMAX ≈ 1000.0 * 2π / 0.0976
        @test b2.TMAX ≈ 2 * b1.TMAX rtol = 1e-9

        # Fixed TMAX ignores ω_ref entirely.
        fixed_spec = parse_spec(write_toml(dir, "ft.toml"; simulation = Dict("TMAX" => 500.0)))
        @test resolve_budget(fixed_spec, 64, 0.0976).TMAX == 500.0
    end
end

@testset "n_samples keeps the sample count constant across N" begin
    mktempdir() do dir
        # modal_E grows as N·nt. Holding nt fixed makes it ∝ N instead of ∝ N².
        spec = parse_spec(write_toml(dir, "ns.toml";
                                     simulation = Dict("scaled_t_max" => 1000.0,
                                                       "n_samples" => 500)))
        for (N, ω) in ((64, 0.0976), (256, 0.0244))
            b  = resolve_budget(spec, N, ω)
            nt = b.TMAX / (spec.DT * b.save_every)
            @test isapprox(nt, 500; rtol = 0.02)
        end
    end
end

@testset "resolve_budget rejects a missing time budget and a bad ω_ref" begin
    mktempdir() do dir
        spec = parse_spec(write_toml(dir, "ok.toml"))
        @test_throws ErrorException resolve_budget(spec, 64, 0.0)
        @test_throws ErrorException resolve_budget(spec, 64, -1.0)

        no_time = parse_spec(write_toml(dir, "nt.toml"; simulation = Dict("TMAX" => 0.0)))
        @test_throws ErrorException resolve_budget(no_time, 64, 0.0976)
    end
end

# ── CaseSetup ─────────────────────────────────────────────────────────────────

@testset "build_case wires system, modes and budget together" begin
    mktempdir() do dir
        spec = parse_spec(write_toml(dir, "c.toml";
                          physics = Dict("N" => 32, "system_type" => "springs",
                                         "nonlinear" => "alpha")))
        case = build_case(spec, 32, 0.1, 0.2)

        # springs ⇒ delta goes to delta_k; alpha ⇒ param goes to alpha.
        @test case.sp.delta_k == 0.2
        @test case.sp.delta_m == 0.0
        @test case.sp.alpha   == 0.1
        @test case.sp.beta    == 0.0

        # Under "mode" both roles coincide, and PBC skips the translation.
        @test case.excited  == [2]
        @test case.ref_mode == 2
        @test case.omega_ref == case.freq[2]
        @test case.budget.TMAX == 100.0
    end
end

@testset "masses/beta route delta and param to the other fields" begin
    mktempdir() do dir
        spec = parse_spec(write_toml(dir, "mb.toml";
                          physics = Dict("system_type" => "masses", "nonlinear" => "beta")))
        case = build_case(spec, 32, 0.3, 0.4)
        @test case.sp.delta_m == 0.4
        @test case.sp.delta_k == 0.0
        @test case.sp.beta    == 0.3
        @test case.sp.alpha   == 0.0
    end
end

@testset "band_ensemble excites a band and clocks off its lowest mode" begin
    mktempdir() do dir
        spec = parse_spec(write_toml(dir, "be.toml";
                          physics = Dict("N" => 32, "init_type" => "band_ensemble",
                                         "branch" => "acoustic", "n_real" => 2)))
        case = build_case(spec, 32, 0.1, 0.1)

        # The whole acoustic band is excited (2:N/2), not a single mode.
        @test length(case.excited) > 1
        @test case.excited == collect(2:16)
        # ref_mode is the band's lowest mode — the two roles genuinely differ here.
        @test case.ref_mode == 2
        @test case.omega_ref == case.freq[2]
        # Bands partition into acoustic and optical for interband diagnostics.
        @test case.k_ac == collect(2:16)
        @test case.k_opt == collect(17:32)
    end
end

@testset "init_mode under band_ensemble is reported, not silently dropped" begin
    mktempdir() do dir
        spec = parse_spec(write_toml(dir, "bi.toml";
                          physics = Dict("init_type" => "band_ensemble",
                                         "branch" => "acoustic", "init_mode" => 7)))
        case = @test_logs (:info,) match_mode = :any build_case(spec, 32, 0.1, 0.1)
        # The band still decides; the point is that it says so.
        @test case.ref_mode == 2
    end
end

@testset "initial_condition puts the requested energy in the chain" begin
    mktempdir() do dir
        spec = parse_spec(write_toml(dir, "ic.toml";
                          physics = Dict("N" => 32, "initial_energy" => 0.5)))
        case   = build_case(spec, 32, 0.1, 0.1)
        q0, v0 = initial_condition(case, spec, spec.seed_base)

        @test length(q0) == 32
        @test v0 == zeros(32)

        # Mode init is deterministic: the seed must not matter.
        q1, _ = initial_condition(case, spec, spec.seed_base + 999)
        @test q1 == q0

        # Harmonic energy of the excited mode: ½ω²A² with A the modal amplitude.
        U  = Diagonal(1.0 ./ sqrt.(case.m)) * case.V
        A  = (U \ q0)[case.ref_mode]
        @test 0.5 * case.freq[case.ref_mode]^2 * A^2 ≈ 0.5 rtol = 1e-8
    end
end

@testset "band_ensemble initial conditions are seed-reproducible" begin
    mktempdir() do dir
        spec = parse_spec(write_toml(dir, "bs.toml";
                          physics = Dict("N" => 32, "init_type" => "band_ensemble",
                                         "branch" => "acoustic")))
        case = build_case(spec, 32, 0.1, 0.1)

        a1, b1 = initial_condition(case, spec, 7)
        a2, b2 = initial_condition(case, spec, 7)
        a3, _  = initial_condition(case, spec, 8)

        @test a1 == a2 && b1 == b2      # same seed ⇒ bit-identical
        @test a1 != a3                  # different seed ⇒ different draw
    end
end

@testset "a config asking for the PBC zero mode fails loudly at build_case" begin
    mktempdir() do dir
        # This is the end-to-end version of the guard: a TOML can request mode 1
        # under PBC, and the failure must surface here rather than as a silently
        # collapsed time axis several CPU-hours later.
        spec = parse_spec(write_toml(dir, "bad.toml";
                          physics = Dict("boundary" => "periodic", "init_mode" => 1)))
        @test_throws ErrorException build_case(spec, 32, 0.1, 0.1)
    end
end

# ── SweepDriver ───────────────────────────────────────────────────────────────

@testset "sweep_tasks builds the full grid, heaviest first" begin
    mktempdir() do dir
        spec = parse_spec(write_toml(dir, "g.toml";
                          physics = Dict("N_values" => [32, 128, 64],
                                         "param_values" => [0.1, 0.2],
                                         "delta_values" => [0.05, 0.1, 0.3])))
        tasks = sweep_tasks(spec)

        @test length(tasks) == 3 * 2 * 3
        @test length(unique(tasks)) == length(tasks)
        # Descending N so a dynamic schedule starts the expensive cases first.
        @test issorted([t.N for t in tasks], rev = true)
        @test Set(t.N for t in tasks) == Set([32, 64, 128])
    end
end

@testset "run_sweep saves results and survives a failing task" begin
    mktempdir() do dir
        out = joinpath(dir, "out.jld2")
        tasks = [1, 2, 3, 4]

        # Task 3 throws; the sweep must keep the other three rather than lose all.
        path = run_sweep(tasks; outfile = out, config_path = "cfg.toml",
                         parallel = :serial) do t
            t == 3 && error("boom")
            (; value = t^2)
        end

        @test path == out
        @test isfile(out)
        saved = load(out)
        @test saved["config"] == "cfg.toml"
        @test length(saved["results"]) == 3
        @test sort([r.value for r in saved["results"]]) == [1, 4, 16]
    end
end

@testset "run_sweep parallel modes agree and bad modes are rejected" begin
    mktempdir() do dir
        expected = [i^2 for i in 1:8]
        for mode in (:serial, :threads, :dynamic)
            out = joinpath(dir, "p_$mode.jld2")
            run_sweep(1:8; outfile = out, config_path = "c.toml", parallel = mode) do t
                (; value = t^2)
            end
            @test sort([r.value for r in load(out)["results"]]) == expected
        end

        @test_throws ErrorException run_sweep(identity, 1:3;
            outfile = joinpath(dir, "x.jld2"), config_path = "c.toml", parallel = :nope)
        @test_throws ErrorException run_sweep(identity, Int[];
            outfile = joinpath(dir, "y.jld2"), config_path = "c.toml")
    end
end

end # testset Compute pipeline
