"""
    plot_threshold_transfer.jl

Open vs closed acoustic→optical channel, as two panels answering the two
robustness axes a referee asks about:

  (a) system size      — x = N,   at the paper's energy density (α = 0.1)
  (b) energy density   — x = λ̄² = α²ε, at N = 64

Each panel has two curves: Δκ = 0.1 (open channel, below the threshold) and
Δκ = 0.7 (closed, above it). The ordinate is E_opt/E, the fraction of energy in
the optical band, read at a fixed `--at` in scaled time.

Why two Δκ instead of a Δκ sweep: with a single-mode initial condition on a
finite lattice, the fast three-wave channel is open only where the discrete
k-grid offers a near-resonant triad with the excited mode, which happens at
Δκ = 0.1 and not at 0.2–0.4. A sweep therefore shows scatter, not an edge at
1/2; the contrast that is robust is open vs closed.

Δκ = 1/2 is the *marginal* case (2ω_ac(π) = ω_op(0) exactly, and k = 0, π are
grid points for every even N), so it is left out; Δκ ≥ 0.6 is genuinely isolated.

Usage:
  julia --project=. scripts/plot/plot_threshold_transfer.jl [--at T] [--open D] [--closed D] <dir_or_file>...
"""

using JLD2, LaTeXStrings, Printf
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils
include("../../src/plot_style.jl");     using .PlotStyle
using Plots

const USAGE = "Usage: julia --project=. scripts/plot/plot_threshold_transfer.jl [--at T] [--open D] [--closed D] <dir_or_file>..."

"""Every `.jld2` under `path`, or `path` itself when it is already one."""
expand(path) = isdir(path) ? [joinpath(path, f) for f in readdir(path) if endswith(f, ".jld2")] :
                             [path]

function parse_args(args)
    at, open_d, closed_d, paths = 1e5, 0.1, 0.7, String[]
    i = 1
    while i <= length(args)
        if     args[i] == "--at";     at       = parse(Float64, args[i+1]); i += 2
        elseif args[i] == "--open";   open_d   = parse(Float64, args[i+1]); i += 2
        elseif args[i] == "--closed"; closed_d = parse(Float64, args[i+1]); i += 2
        else push!(paths, args[i]); i += 1
        end
    end
    isempty(paths) && error(USAGE)
    at, open_d, closed_d, paths
end

"""
    optical_fraction(res, at) -> Float64

E_opt/E at scaled time `at`, or `NaN` when the run ends before it. Under
`init_type="mode"` the optical band is the upper half of the spectrum
(`CaseSetup.bands`), which is what `k_opt` holds for those runs.
"""
function optical_fraction(res, at::Float64)
    t = Float64.(Vector(res.scaled_t))
    t[end] < at * 0.98 && return NaN   # 2% slack: n_samples rarely lands on `at`
    E = Float64.(Matrix(res.modal_E))
    size(E, 1) > size(E, 2) && (E = E')
    N = size(E, 1)
    i = min(searchsortedfirst(t, at), length(t))
    sum(@view E[N÷2+1:end, i]) / sum(@view E[:, i])
end

"""
All results in `paths`, as `(N, param, Delta, frac)` rows at fixed `at`. A
(N, α, Δκ) case seen twice keeps its first row, so overlapping sweeps (a short
run and its longer rerun) don't draw duplicate points.
"""
function rows(paths, at)
    out = NamedTuple[]
    seen = Set{Tuple{Int,Float64,Float64}}()
    for p in paths, f in expand(p), r in load(f)["results"]
        key = (r.N, Float64(r.param), Float64(r.Delta))
        key in seen && continue
        push!(seen, key)
        push!(out, (N = r.N, param = key[2], Delta = key[3], frac = optical_fraction(r, at)))
    end
    out
end

const EPS = 0.00695   # energy density of every run in these sweeps

"""Open/closed pair of curves: x from `xof(row)`, one point per row at that Δκ."""
function pair(rs, xof, open_d, closed_d)
    cs = Curve[]
    for (j, (d, name)) in enumerate(((open_d, "open"), (closed_d, "closed")))
        sel = sort(filter(r -> r.Delta == d && isfinite(r.frac), rs), by = xof)
        isempty(sel) && continue
        push!(cs, Curve([xof(r) for r in sel], [r.frac for r in sel],
                        latexstring("\\Delta\\kappa = $d\\ (\\mathrm{$name})");
                        color = j == 1 ? "#1a3a6b" : "#d62728",
                        linestyle = j == 1 ? :solid : :dash, linewidth = 2.6))
        @printf("  Δκ=%-4.2f %s\n", d,
                join((@sprintf("%g:%.1e", xof(r), r.frac) for r in sel), "  "))
    end
    cs
end

function panel(cs, xlabel, title, xticks_)
    p = plot(; xlabel = xlabel, ylabel = L"E_{\mathrm{opt}}/E", title = title,
             xscale = :log10, yscale = :log10, xticks = xticks_,
             framestyle = :box, grid = false, legend = :topleft, ylims = (1e-5, 1.5))
    draw!(p, cs)
    for c in cs
        scatter!(p, c.x, c.y; color = c.color, markerstrokewidth = 0, ms = 6, label = "")
    end
    p
end

function main()
    at, open_d, closed_d, paths = parse_args(ARGS)
    apply_style!()
    rs = rows(paths, at)
    isempty(rs) && error("No results found in: " * join(paths, ", "))

    println("Panel (a) — size, α = 0.1:")
    a = panel(pair(filter(r -> r.param == 0.1, rs), r -> Float64(r.N), open_d, closed_d),
              L"N", L"(a)\ \alpha = 0.1",
              ([64, 128, 256], ["64", "128", "256"]))

    println("Panel (b) — energy density, N = 64:")
    lam2(r) = r.param^2 * EPS
    sel = filter(r -> r.N == 64 && r.param <= 0.2 + 1e-9, rs)   # α = 0.3 is outside the weak-coupling regime
    xs = sort(unique(lam2.(sel)))
    b = panel(pair(sel, lam2, open_d, closed_d),
              L"\bar{\lambda}^{2} = \alpha^{2}\epsilon", L"(b)\ N = 64",
              (xs, [@sprintf("%.1f", x * 1e5) for x in xs]))
    xlabel!(b, L"\bar{\lambda}^{2} = \alpha^{2}\epsilon\ (\times 10^{-5})")

    fig = plot(a, b; layout = (1, 2), size = (1200, 520),
               left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
    save_fig(fig, "results/figures/threshold",
             @sprintf("threshold_open_closed_at%.0e", at))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
