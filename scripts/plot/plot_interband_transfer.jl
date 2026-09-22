"""
    plot_interband_transfer.jl

Interband (acoustic→optical) energy transfer for an open channel (Δκ below the
1/2 threshold) against a closed one (Δκ above it), as time series rather than
as a value read at one instant:

  top row     — system size,    N = 64, 128, 256 at the paper's energy density
  bottom row  — energy density, λ̄² = α²ε at N = 64

One small panel per (N or λ̄²), each holding a single open/closed pair in the
same two colours throughout, so no panel ever asks the reader to tell three
shades apart; the widening gap from left to right *is* the result. The guide at
E_opt/E = 10⁻² marks appreciable transfer: no closed-channel run reaches it
within 10⁶ cycles, at any size or density explored.

Why E_opt/E and not the spectral entropy used elsewhere in the repo: S̄ counts
*how many* modes carry energy, and the acoustic band alone is half the
spectrum, so intra-acoustic spreading alone drives S̄/ln N to ≈0.5 without a
single quantum crossing the gap. Measured at α=0.1, the open/closed separation
is 1.02× in S̄ at N=128 but 25× in E_opt/E. The threshold gates *how much*
energy crosses, which is what this observable tracks.

The initial plateau at ~1e-4 is not transfer: it is the non-resonant distortion
of the single-mode initial condition projected onto the optical modes, and it
scales as λ̄². Transfer is the sustained growth above it.

Usage:
  julia --project=. scripts/plot/plot_interband_transfer.jl [--open D] [--closed D] <dir_or_file>...
"""

using JLD2, LaTeXStrings, Printf
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils
include("../../src/plot_style.jl");     using .PlotStyle
using Plots

const USAGE = "Usage: julia --project=. scripts/plot/plot_interband_transfer.jl [--open D] [--closed D] <dir_or_file>..."
const EPS = 0.00695            # energy density shared by every run plotted here
const ALPHA_PAPER = 0.1        # the manuscript's nonlinearity
const N_PANEL_B = 64
# α = 0.05 is kept as the leftmost density panel precisely because nothing happens
# there within 1e6 cycles: a flat pair next to a rising one shows that the window
# is long enough to resolve transfer when there is any. α = 0.15 is dropped as an
# intermediate step between 0.1 and 0.2; α = 0.3 is left out because at ε_eff ≈
# 0.06 the system is no longer weakly coupled and the threshold blurs.
const ALPHAS_B = (0.05, 0.1, 0.2)
const GUIDE = 1e-2             # "appreciable transfer" level
const C_OPEN   = "#2a78d6"     # one hue pair for every panel: identity comes from
const C_CLOSED = "#eb6834"     # the panel title, never from the colour
# Linestyle repeats the channel distinction: in greyscale these two hues sit at
# 119 and 143 (1.38:1), so colour alone does not survive a black-and-white print.
const S_OPEN, S_CLOSED = :solid, :dash
# Weight repeats it once more: a legend swatch is too short to show a dash
# pattern, but thick-vs-thin reads at any size and in black and white.
const W_OPEN, W_CLOSED = 2.8, 1.2

"""λ̄² as `m×10^e` with the exponent taken from the value, so a panel at 10⁻⁵ is
not labelled in units of 10⁻⁴."""
function sci_label(x::Real)
    e = floor(Int, log10(x))
    latexstring(@sprintf("\\bar{\\lambda}^2 = %.1f\\times10^{%d}", x / 10.0^e, e))
end

expand(p) = isdir(p) ? [joinpath(p, f) for f in readdir(p) if endswith(f, ".jld2")] : [p]

function parse_args(args)
    open_d, closed_d, paths = 0.1, 0.7, String[]
    i = 1
    while i <= length(args)
        if     args[i] == "--open";   open_d   = parse(Float64, args[i+1]); i += 2
        elseif args[i] == "--closed"; closed_d = parse(Float64, args[i+1]); i += 2
        else push!(paths, args[i]); i += 1
        end
    end
    isempty(paths) && error(USAGE)
    open_d, closed_d, paths
end

"""
    series(res) -> (t, f)

Optical-band energy fraction against scaled time, thinned to one point per
logarithmic bin. Under `init_type="mode"` the optical band is the upper half of
the spectrum, matching `CaseSetup.bands`.
"""
function series(res)
    t = Float64.(Vector(res.scaled_t))
    E = Float64.(Matrix(res.modal_E))
    size(E, 1) > size(E, 2) && (E = E')
    N = size(E, 1)
    keep = findall(>(1.0), t)                     # log x: drop t ≤ 1
    f = vec(sum(@view(E[N÷2+1:end, keep]), dims = 1)) ./ vec(sum(@view(E[:, keep]), dims = 1))
    logdownsample(t[keep], f)
end

"""Every run in `paths` reaching 1e6 cycles, keyed by `(N, α, Δκ)`, first wins."""
function load_runs(paths)
    runs = Dict{Tuple{Int,Float64,Float64},Any}()
    for p in paths, file in expand(p), r in load(file)["results"]
        key = (r.N, Float64(r.param), Float64(r.Delta))
        haskey(runs, key) && continue
        last(r.scaled_t) < 9e5 && continue        # a shorter window would end mid-figure
        runs[key] = r
    end
    runs
end

"""Open/closed curves for one panel; only the first panel shows the key."""
function pair(runs, key, open_d, closed_d)
    cs = Curve[]
    for (d, color, style, w) in ((open_d, C_OPEN, S_OPEN, W_OPEN),
                                 (closed_d, C_CLOSED, S_CLOSED, W_CLOSED))
        r = get(runs, key(d), nothing)
        isnothing(r) && (@printf("  ⚠ missing: %s Δκ=%.2f\n", string(key(d)), d); continue)
        t, f = series(r)
        push!(cs, Curve(t, f, latexstring("\\Delta\\kappa = $d");
                        color = color, linestyle = style, linewidth = w))
        @printf("  %-22s Δκ=%.2f  end=%.2e\n", string(key(d)), d, f[end])
    end
    cs
end

"""One small panel: axes only on the outer edges, so the 2×3 grid reads as a unit."""
function panel(cs, title; ylab::Bool, xlab::Bool, legend::Bool = false)
    p = plot(; title = title, titlefontsize = 12,
             xlabel = xlab ? L"t" : "", ylabel = ylab ? L"E_{\mathrm{opt}}/E" : "",
             xscale = :log10, yscale = :log10,
             xticks = (10.0 .^ (2:2:6), [latexstring("10^{$i}") for i in 2:2:6]),
             yticks = (10.0 .^ (-6:2:0), [latexstring("10^{$i}") for i in -6:2:0]),
             xlims = (10.0, 1.2e6), ylims = (1e-7, 1.0),
             framestyle = :box, grid = false,
             legend = legend ? :topleft : false, legendfontsize = 10,
             top_margin = 5Plots.mm, bottom_margin = 3Plots.mm)
    draw!(p, cs)
    hline!(p, [GUIDE]; color = GRAY_GUIDE, linestyle = :dot, lw = 1.2, label = "")
    p
end

function main()
    open_d, closed_d, paths = parse_args(ARGS)
    apply_style!()
    runs = load_runs(paths)
    isempty(runs) && error("No 1e6-cycle runs found in: " * join(paths, ", "))

    panels = Any[]
    println("Top row — size, α = $ALPHA_PAPER:")
    for (j, N) in enumerate((64, 128, 256))
        cs = pair(runs, d -> (N, ALPHA_PAPER, d), open_d, closed_d)
        push!(panels, panel(cs, latexstring("N = $N");
                            ylab = j == 1, xlab = false, legend = j == 1))
    end
    println("Bottom row — energy density, N = $N_PANEL_B:")
    for (j, a) in enumerate(ALPHAS_B)
        cs = pair(runs, d -> (N_PANEL_B, a, d), open_d, closed_d)
        push!(panels, panel(cs, sci_label(a^2 * EPS); ylab = j == 1, xlab = true))
    end

    fig = plot(panels...; layout = (2, 3), size = (1250, 760),
               left_margin = 6Plots.mm, bottom_margin = 8Plots.mm)
    save_fig(fig, "results/figures/threshold", "interband_transfer")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
