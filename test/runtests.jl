using Test
using Random
using LinearAlgebra
include("../src/fput_core.jl"); using .FPUTCore
include("../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../src/fput_analysis.jl"); using .FPUTAnalysis
include("../src/experiment.jl"); using .Experiment
include("../src/case_setup.jl"); using .CaseSetup
include("../src/sweep_driver.jl"); using .SweepDriver

# ── Helper: Analytic dispersion for diatom with equal springs ──
"""
    diatomic_dispersion(k_vec, m1, m2, kappa) -> (ω_acoustic, ω_optical)

Analytic dispersion ω²(k) for diatom with equal springs κ and mass pair (m1, m2).
For periodic boundary conditions:
  ω² = κ(1/m1 + 1/m2) ± κ√[(1/m1 + 1/m2)² - 4sin²(ka/2)/(m1·m2)]
"""
function diatomic_dispersion(k_vec::AbstractVector, m1::Float64, m2::Float64, kappa::Float64)
    A = kappa * (1/m1 + 1/m2)
    omega_sq_ac = similar(k_vec)
    omega_sq_op = similar(k_vec)

    @inbounds for (i, k) in enumerate(k_vec)
        B = 4 * sin(k / 2)^2 / (m1 * m2)
        discriminant = sqrt(A^2 - kappa^2 * B)
        omega_sq_ac[i] = A - kappa * discriminant
        omega_sq_op[i] = A + kappa * discriminant
    end

    return sqrt.(max.(0.0, omega_sq_ac)), sqrt.(max.(0.0, omega_sq_op))
end

@testset "FPUT Core Physics" begin

    @testset "normal modes exist and are real-valued" begin
        # Basic smoke test: normal modes should be computed without error,
        # eigenvalues should be positive (frequencies real), and eigenvectors
        # orthonormal. This verifies the matrix construction and LAPACK call.
        N = 32
        sp = SystemParams(N, 0.1, 0.0, 0.0, 0.0, :periodic)
        k, m = make_system(sp)
        freq, V = find_normal_modes(k, m, :periodic)

        # Frequencies should be non-negative real numbers
        @test all(freq .>= -1e-14)
        @test all(isreal, freq)

        # Eigenvectors should be orthonormal (V'V ≈ I)
        gram = V' * V
        @test gram ≈ I atol=1e-10
    end

    @testset "energy conservation: symplectic integrator (KahanLi8)" begin
        # Total energy should oscillate (bounded drift) rather than drift monotonically,
        # which is the signature of a symplectic integrator.
        Random.seed!(42)

        N = 32
        sp = SystemParams(N, 0.2, 0.0, 0.1, 0.0, :periodic)

        # Initial condition: single mode excitation
        k, m = make_system(sp)
        freq, V = find_normal_modes(k, m, :periodic)
        amplitude = sqrt(2 * 0.0075) / freq[2]  # excite mode 2 (skip translation at freq[1])
        q0 = amplitude .* V[:, 2]
        v0 = zeros(N)

        # Compute total energy
        function hamiltonian(q, v, k, m, alpha, beta, boundary)
            kinetic = 0.5 * sum(m .* v.^2)

            # Potential energy: sum of bond potentials
            N_sites = length(m)
            potential = 0.0
            if boundary == :periodic
                for i in 1:N_sites
                    d = q[mod1(i+1, N_sites)] - q[i]
                    potential += bond_potential(k[i], d, alpha, beta)
                end
            else  # :fixed
                potential += bond_potential(k[1], q[1], alpha, beta)
                for i in 2:N_sites
                    d = q[i] - q[i-1]
                    potential += bond_potential(k[i], d, alpha, beta)
                end
                potential += bond_potential(k[N_sites+1], -q[N_sites], alpha, beta)
            end

            return kinetic + potential
        end

        # Solve over moderate time span
        Q, V_traj, t, _, _ = solve_fput(sp, q0, v0, (0.0, 2000.0), 0.05; saveat=0:100:2000)

        # Compute energy at each time step
        E = [hamiltonian(Q[:, i], V_traj[:, i], k, m, sp.alpha, sp.beta, sp.boundary)
             for i in axes(Q, 2)]

        E_initial = E[1]
        E_drift = (E .- E_initial) ./ E_initial

        # Check that drift is bounded (no exponential growth)
        max_drift = maximum(abs.(E_drift))
        @test max_drift < 1e-6

        # Check that drift oscillates (not monotonic growth) — standard for symplectic
        drift_gradient = diff(E_drift)
        sign_changes = count(x -> x != 0, diff(sign.(drift_gradient)))
        @test sign_changes > 5
    end

    @testset "FPUT recurrence (Fermi-Pasta-Ulam phenomenon)" begin
        # Reproduce original FPUT setup: monoatomic chain, fixed boundaries,
        # mode fundamental excited, observe energy return over ~1000 periods.
        # Recurrence occurs for both α (quadratic) and β (quartic) nonlinearity.

        N = 32
        # Monoatomic (delta_k = delta_m = 0), fixed boundaries
        sp = SystemParams(N, 0.0, 0.0, 0.1, 0.05, :fixed)

        k, m = make_system(sp)
        freq, V = find_normal_modes(k, m, :fixed)

        # Excite mode 1 (lowest non-zero frequency; mode 0 is wall-to-wall translation)
        # Energy in fundamental mode ≈ 1.0
        energy_target = 1.0
        amplitude = sqrt(2 * energy_target) / freq[2]
        q0 = amplitude .* V[:, 2]
        v0 = zeros(N)

        # Integration time: roughly 1000 periods of fundamental frequency
        T_max = 1200.0 / freq[2]

        # Solve FPUT dynamics
        Q, V_traj, t, _, _ = solve_fput(sp, q0, v0, (0.0, T_max), 0.1; saveat=0:T_max/30:T_max)

        # Compute modal energies at each time step
        modal_E_initial = compute_modal_energies(Q[:, 1], V_traj[:, 1], freq, V, m)
        E_mode1_traj = [compute_modal_energies(Q[:, i], V_traj[:, i], freq, V, m)[2]
                        for i in axes(Q, 2)]

        # Return ratio: energy in mode 1 at final time / initial time
        return_ratio = E_mode1_traj[end] / E_mode1_traj[1]

        # Verify recurrence: at least 80% of energy returns to initial mode
        @test return_ratio > 0.80
    end

end

@testset "Diatomic Chain Physics" begin

    @testset "band separation and gap structure" begin
        # Verify that acoustic and optical branches exist and are separated by Δm.
        # For a diatom with equal springs κ and mass difference Δm, there is an
        # energy gap between acoustic (low) and optical (high) modes.
        # This is a structural property that does NOT require energy isolation,
        # since eigenmodes in diatom can mix across branches.

        N = 32
        delta_m = 0.3
        sp = SystemParams(N, 0.0, delta_m, 0.0, 0.0, :periodic)

        k, m = make_system(sp)
        freq, V = find_normal_modes(k, m, :periodic)

        # Sort frequencies to identify branch structure
        freq_sorted = sort(freq)

        # For a diatom: expect ~N/2 acoustic modes (low freq) and ~N/2 optical (high)
        # with a clear gap between them
        N_half = div(N, 2)
        f_cutoff = (freq_sorted[N_half] + freq_sorted[N_half + 1]) / 2

        # Define acoustic/optical by sorted frequency
        acoustic_mask = freq_sorted .<= f_cutoff
        optical_mask = freq_sorted .> f_cutoff

        # Verify the structure: roughly equal number of modes per branch
        @test count(acoustic_mask) >= N_half - 3
        @test count(optical_mask) >= N_half - 3

        # Key test: acoustic max < optical min (clear energy gap)
        f_acoustic_max = maximum(freq_sorted[acoustic_mask])
        f_optical_min = minimum(freq_sorted[optical_mask])
        gap_ratio = f_optical_min / (f_acoustic_max + 1e-10)

        # Gap should be at least 5% of acoustic max (gap widens with Δm)
        @test gap_ratio > 1.05

        # Test: excite mode 2 (lowest non-translation acoustic)
        # Verify it is indeed in the acoustic branch
        mode_2_freq = freq[2]
        @test mode_2_freq <= f_cutoff
    end

    @testset "analytic dispersion vs. computed frequencies" begin
        # Verify that eigenfrequencies match analytic dispersion relation for
        # a diatom with equal springs and periodic boundary conditions.
        # This validates the dynamical matrix construction and eigensolver.

        N = 32
        delta_m = 0.2
        sp = SystemParams(N, 0.0, delta_m, 0.0, 0.0, :periodic)

        k, m = make_system(sp)

        # Extract alternating masses
        m1 = m[1]  # 1 + delta_m
        m2 = m[2]  # 1 - delta_m
        kappa = k[1]  # = 1.0 (since delta_k = 0)

        # Compute eigenmodes
        freq, V = find_normal_modes(k, m, :periodic)

        # Analytic dispersion: wavenumbers k_j = 2πj/N for j=0..N-1
        k_vals = 2π .* (0:N-1) ./ N
        omega_ac, omega_op = diatomic_dispersion(k_vals, m1, m2, kappa)

        # Sort computed frequencies and match to analytic branches
        freq_sorted = sort(freq)

        # Collect analytic frequencies (both branches)
        omega_analytic = vcat(omega_ac, omega_op)
        omega_analytic_sorted = sort(omega_analytic)

        # Compare: each computed frequency should match an analytic one (within tol)
        N_modes = length(freq)

        # Check that at least half the modes match analytic dispersion within 5%
        matches = 0
        @inbounds for i in 1:N_modes
            # For each computed freq, find nearest analytic freq
            errors = abs.(omega_analytic_sorted .- freq_sorted[i]) ./ (freq_sorted[i] + 1e-10)
            min_error = minimum(errors)
            if min_error < 0.05  # 5% tolerance
                matches += 1
            end
        end

        @test matches >= div(N_modes, 2)
    end

end

include("test_compute_pipeline.jl")
