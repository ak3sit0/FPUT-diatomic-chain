using Test
using Random
using LinearAlgebra
include("../src/fput_core.jl"); using .FPUTCore
include("../src/fput_fast_runner.jl"); using .FPUTFastRunner
include("../src/fput_analysis.jl"); using .FPUTAnalysis

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

end
