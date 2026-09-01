"""
    plot_ftmle.jl

Single-panel figure: λ(t) in log-log, one curve per input file, plus a fitted
power law t^{-δ} over the reference case (Δκ = `DELTA_REF`), fit on the central
50% of that curve.

Usage:
  julia --project=. scripts/plot/plot_ftmle.jl <ftmle1.jld2> [ftmle2.jld2 ...]
"""

using JLD2, LaTeXStrings, Statistics, Dates
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils
include("../../src/plot_style.jl");     using .PlotStyle

const USAGE = "Usage: julia --project=. scripts/plot/plot_ftmle.jl <ftmle1.jld2> [ftmle2.jld2 ...]"
const DELTA_REF = 0.05   # case the power-law reference is fit on

"""
    fit_powerlaw(t, lam; t_max=1e4, n_bins=20) -> (c, δ)

Log-binned linear regression in log-log: split `[t[1], t_max]` into `n_bins`
equal-width bins in `log(t)`, average λ within each, fit on the bin means. Every
decade weighs equally regardless of point density.
"""
function fit_powerlaw(t, lam; t_max=1e4, n_bins=20)
    # Log-space binning ensures each decade contributes equally to the fit, avoiding
    # the usual problem where later (sparser) time points dominate linear regression.
    mask = (t .<= t_max) .& (t .> 0) .& (lam .> 0)
    t_m, l_m = t[mask], lam[mask]
    isempty(t_m) && error("No data in the fit range")

    # Bin edges in log space (equal width in log-scale, hence many time points per bin at early times).
    edges = exp10.(range(log10(t_m[1]), log10(t_max); length=n_bins+1))
    t_bin, l_bin = Float64[], Float64[]
    for i in 1:n_bins
        idx = findall(s -> edges[i] <= s < edges[i+1], t_m)
        isempty(idx) && continue
        # Geometric mean (mean of logs) for both axes: better center-of-mass in log-log space.
        push!(t_bin, exp(mean(log.(t_m[idx]))))
        push!(l_bin, exp(mean(log.(l_m[idx]))))
    end

    # Linear regression in log-log space: log(λ) = log(c) + (-δ)·log(t).
    x, y = log.(t_bin), log.(l_bin)
    mx, my = mean(x), mean(y)
    slope = sum((x .- mx) .* (y .- my)) / sum((x .- mx).^2)  # standard least-squares slope
    exp(my - slope * mx), -slope  # (c, δ) such that λ ≈ c·t^{-δ}
end

"""One λ(t) curve per file, plus the c·t^{-δ} reference fit on `DELTA_REF`."""
function ftmle_curves(paths)
    curves = Curve[]
    ref = nothing
    for (k, path) in enumerate(paths)
        d = load(path)
        t, lam, delta = Float64.(d["t_cycles"]), Float64.(d["lambda"]), Float64(d["delta_k"])
        # Filter to positive, finite values (NaN and ≤0 cannot be plotted on log axes).
        valid = findall(x -> x > 0 && isfinite(x), lam)
        t, lam = t[valid], lam[valid]

        is_ref = isapprox(delta, DELTA_REF; atol=1e-10)  # store fit only for the reference Δκ
        is_ref && (ref = (t[1], t[end], fit_powerlaw(t, lam)...))

        t_ds, l_ds = logdownsample(t, lam)
        push!(curves, Curve(t_ds, l_ds, latexstring("\\Delta\\kappa = $delta");
                            cyclic(k; palette = FTMLE_PALETTE, styles = FTMLE_LINESTYLES,
                                   lw = is_ref ? 3.5 : 1.5)...))
        println("  ✓ Δκ=$delta ($(length(t)) points)")
    end

    if !isnothing(ref)
        t0, t1, c, δ = ref
        t_fit = exp10.(range(log10(t0), log10(t1); length=400))
        push!(curves, Curve(t_fit, c .* t_fit .^ (-δ),
                            latexstring("c\\,t^{-$(round(δ; digits=3))}");
                            color = :black, linestyle = :solid, linewidth = 2.0))
        println("Fit (Δκ=$DELTA_REF): c = $(round(c; sigdigits=3))  δ = $(round(δ; digits=4))")
    end
    curves
end

function boundary_tag(paths)
    tags = unique(get(load(p), "boundary", "periodic") for p in paths)
    length(tags) == 1 ? tags[1] : "mixed"
end

function main()
    isempty(ARGS) && error(USAGE)
    apply_style!()

    fig = logplot(xlabel = L"t\ \mathrm{(cycles)}", ylabel = L"\lambda(t)",
                  yscale = :log10, legend = :topright,
                  grid = true, gridalpha = 0.2, gridstyle = :dot,
                  size = (800, 500))
    draw!(fig, ftmle_curves(ARGS))

    deltas = sort(unique(Float64(load(p)["delta_k"]) for p in ARGS))
    tag = join(replace.(string.(deltas), "." => "p"), "_")
    save_fig(fig, "results/figures/liapunov_exponent",
             "ftmle_$(boundary_tag(ARGS))_delta$(tag)_$(today())"; exts = (".pdf",))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
