"""
    compute_scattering_rate.jl

Compute the effective acoustic-to-optical scattering rate

    R(η) = ∫ dk₁ dk₂ |Γ_aao(k₁,k₂)|² δ[ω₋(k₁)+ω₋(k₂)−ω₊(k₃)]

as a line integral over the resonance manifold C_η, using the co-area formula:

    R(η) = ∫_{C_η} |Γ_aao|² / |∇residual| d𝓁

Requires: Contour.jl, Interpolations.jl, CairoMakie.jl, LaTeXStrings.jl
Include plot_gamma_with_resonance.jl for compute_gamma, resonance_matrix, omega_branch.
"""

using LinearAlgebra, CairoMakie, LaTeXStrings
using Contour: contours, levels, lines as contour_lines, coordinates
using Interpolations: interpolate, BSpline, Cubic, Line, OnGrid, scale

include("../../src/fput_coupling.jl")
using .FPUTCoupling

# ── Dispersion derivatives (analytic) ─────────────────────────────────────────

"""
    domega_dk(k, kA, kB, branch)

Analytic group velocity dω/dk for acoustic (branch=1) or optical (branch=2).
From differentiating ω²(k) = kA+kB ± √((kA+kB)²−4kA·kB·sin²(k/2)).
"""
function domega_dk(k, kA, kB, branch)
    k    = mod(k + π, 2π) - π
    ksum = kA + kB
    kprod = kA * kB
    sin2  = sin(k/2)^2
    disc  = ksum^2 - 4kprod * sin2
    disc  = max(disc, 1e-14)   # avoid sqrt(0) at zone boundary

    # d(disc)/dk = −4kA·kB·sin(k/2)·cos(k/2) = −2kA·kB·sin(k)
    d_disc_dk = -2kprod * sin(k)

    ω = (branch == 1) ? sqrt(max(ksum - sqrt(disc), 0.0)) :
                        sqrt(ksum + sqrt(disc))
    ω < 1e-12 && return 0.0   # ω=0 at Γ point for acoustic branch

    # dω/dk = (1/2ω) · d(ω²)/dk = (1/2ω) · (∓ d(√disc)/dk)
    #        = (1/2ω) · (∓ d_disc_dk / (2√disc))
    sign = (branch == 1) ? -1.0 : +1.0
    return sign * d_disc_dk / (4ω * sqrt(disc))
end

# ── Gradient of the resonance residual ────────────────────────────────────────

"""
    grad_residual(k1, k2, kA, kB)

Gradient of f(k1,k2) = ω₋(k1) + ω₋(k2) − ω₊(k3), with k3 = −k1−k2.
∂f/∂k1 = dω₋/dk|_{k1} + dω₊/dk|_{k3}   (chain rule: ∂k3/∂k1 = −1)
∂f/∂k2 = dω₋/dk|_{k2} + dω₊/dk|_{k3}
"""
function grad_residual(k1, k2, kA, kB)
    k3 = mod(-k1 - k2 + π, 2π) - π
    v1 = domega_dk(k1, kA, kB, 1)
    v2 = domega_dk(k2, kA, kB, 1)
    v3 = domega_dk(k3, kA, kB, 2)   # ∂k3/∂k1 = ∂k3/∂k2 = −1 → sign flip
    return (v1 - v3, v2 - v3)
end

# ── Line integral ──────────────────────────────────────────────────────────────

"""
    scattering_rate(Δκ; Nk=601, Ngrid=401, min_Δκ=1e-3) -> Float64

Evaluate R(Δκ) by:
1. Computing |Γ_aao|² on a (Ngrid×Ngrid) grid.
2. Extracting the resonance curve C via Contour.jl.
3. Integrating |Γ_aao|² / |∇residual| along C by the trapezoid rule.

Returns 0.0 for Δκ < min_Δκ (no physical resonance manifold).
"""
function scattering_rate(Δκ; Nk=601, Ngrid=401, min_Δκ=1e-3)
    Δκ < min_Δκ && return 0.0

    kA = 1.0 + Δκ
    kB = 1.0 - Δκ

    # ── Step 1: |Γ_aao|² on the grid ────────────────────────────────────────
    kplot = range(-π, π, length=Ngrid)
    _, _, Gamma = compute_gamma(kA, kB, 1.0; Nk=Nk, Ngrid=Ngrid)
    # aao is index 2: s1=1(a), s2=1(a), s3=2(o) → (1-1)*4+(1-1)*2+(2-1)+1 = 2
    Gamma2 = abs2.(Gamma[2])   # |Γ_aao|²  — note: abs2, not abs

    # Bilinear interpolant for fast evaluation on arbitrary (k1,k2)
    # Interpolations.jl expects the grid as a range
    itp = interpolate(Gamma2, BSpline(Cubic(Line(OnGrid()))))
    itp_scaled = scale(itp, kplot, kplot)

    # ── Step 2: resonance curve via Contour.jl ───────────────────────────────
    D = resonance_matrix(collect(kplot), kA, kB, 1, 1, 2)
    cs = contours(collect(kplot), collect(kplot), D, [0.0])

    R_total = 0.0

    # Iterate over the contour levels (here we only ask for the [0.0] one)
    for c_level in levels(cs)
        # Extract the lines/arcs of this level using the new name
        for cl in contour_lines(c_level)
            xs, ys = coordinates(cl) # k1, k2 coordinates along the arc
            n = length(xs)
            n < 2 && continue

            # ── Step 3 & 4: trapezoid integration along the arc ─────────────────
            # (The rest of the for i in 1:(n-1) loop stays exactly the same)
            for i in 1:(n-1)
                k1m = 0.5*(xs[i] + xs[i+1])
                k2m = 0.5*(ys[i] + ys[i+1])

                dl = sqrt((xs[i+1]-xs[i])^2 + (ys[i+1]-ys[i])^2)

                # |Γ_aao(k1m, k2m)|² via interpolant
                gamma2_val = itp_scaled(k1m, k2m)

                # |∇residual| at midpoint
                gx, gy = grad_residual(k1m, k2m, kA, kB)
                grad_norm = sqrt(gx^2 + gy^2)
                grad_norm < 1e-12 && continue   # skip degenerate points

                R_total += gamma2_val / grad_norm * dl
            end
        end
    end

    return R_total
end

# ── Sweep over Δκ ─────────────────────────────────────────────────────────────

"""
    sweep_scattering_rate(delta_values; kwargs...) -> Vector{Float64}

Evaluate R(Δκ) for each entry in delta_values.
Prints progress since each call takes a few seconds.
"""
function sweep_scattering_rate(delta_values; Nk=601, Ngrid=401)
    R_vals = zeros(length(delta_values))
    for (i, Δκ) in enumerate(delta_values)
        R_vals[i] = scattering_rate(Δκ; Nk=Nk, Ngrid=Ngrid)
        @show Δκ, R_vals[i]
    end
    return R_vals
end

# ── Plot ──────────────────────────────────────────────────────────────────────

function plot_scattering_rate(delta_values, R_vals; outdir="results/figures/scattering_rate")
    mkpath(outdir)

    # Normalize so the peak = 1 for visual clarity
    R_norm = R_vals ./ max(maximum(R_vals), 1e-30)

    i_max = argmax(R_norm)
    Δκ_opt = delta_values[i_max]

    # Convert Δκ to η = (1 - Δκ)/(1 + Δκ)
    eta_values = @. (1 - delta_values) / (1 + delta_values)
    eta_opt = (1 - Δκ_opt) / (1 + Δκ_opt)

    eta_max = (1 - 0.05) / (1 + 0.05)   # η corresponding to Δκ = 0.05

    fig = Figure(size=(1100, 650))
    ax  = Axis(fig[1, 1],
               #xlabel=L"\eta = \frac{1-\Delta\kappa}{1+\Delta\kappa}",
               xlabel=L"\eta",  # = 1-\Delta\kappa) / 1+\Delta\kappa  ",
               ylabel=L"R(\eta) / R_\mathrm{max}",
               title=L"Effective acoustic-optical scattering rate $R(\eta)$",
               xlabelsize=29, ylabelsize=29, titlesize=32,
               xticklabelsize=22, yticklabelsize=22)

    # Add subtle grid
    hlines!(ax, [0.0, 0.25, 0.5, 0.75, 1.0], color=:gray, alpha=0.2, linewidth=0.5)
    vlines!(ax, collect(0.0:0.1:1.0), color=:gray, alpha=0.2, linewidth=0.5)

    xlims!(ax, 0.0, eta_max)

    # Main curve (dark blue)
    lines!(ax, eta_values, R_norm, linewidth=3.5, color=:darkblue, label=L"R(\eta)")

    # Prominent maximum point (very dark blue)
    scatter!(ax, [eta_opt], [1.0], color=:darkblue, markersize=22, strokewidth=4,
             strokecolor=:white, label=latexstring("\\text{Maximum at } \\eta \\approx $(round(eta_opt, digits=3))"))

    # Vertical line at the maximum (dark blue, dashed)
    vlines!(ax, [eta_opt], linestyle=:dash, color=:darkblue, linewidth=2.0, alpha=0.7)

    # Clean legend
    axislegend(ax, position=:lt, fontsize=16, framevisible=true,
               backgroundcolor=(:white, 0.8), labelsize=24)

    save(joinpath(outdir, "scattering_rate_vs_eta.pdf"), fig)
    return fig
end

# ── Entry point ───────────────────────────────────────────────────────────────
if abspath(PROGRAM_FILE) == @__FILE__
    # Dense sweep from 0.05 to 0.5
    delta_values = collect(0.05:0.01:0.50)

    R_vals = sweep_scattering_rate(delta_values; Nk=601, Ngrid=401)

    # Save raw data
    mkpath("results/scattering_rate")
    open("results/scattering_rate/R_vs_delta.csv", "w") do io
        println(io, "delta_kappa,R")
        for (d, r) in zip(delta_values, R_vals)
            println(io, "$d,$r")
        end
    end

    fig = plot_scattering_rate(delta_values, R_vals)
    display(fig)
end

# Suggested test (add to test/runtests.jl):
# @testset "scattering_rate" begin
#     @test scattering_rate(0.0)   == 0.0          # no manifold at Δκ=0
#     @test scattering_rate(0.51)  == 0.0          # above threshold
#     @test scattering_rate(0.1)    > 0.0          # finite in range
#     @test scattering_rate(0.49)  >= 0.0          # approaching threshold
# end