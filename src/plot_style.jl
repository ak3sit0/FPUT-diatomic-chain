"""
    plot_style.jl

The aesthetic layer for every Plots.jl figure in the repo. A plotting script
should contain no `default(...)`, no `plot(...)` boilerplate, no hand-built
`10^k` ticks and no `savefig` — it builds `Curve`s and calls the four functions
below.

Everything that legitimately differs between figures (`size`, `grid`, `yscale`,
margins) is a keyword of `logplot`, not a fork of it.
"""
module PlotStyle

using Plots, LaTeXStrings
using ..PlottingUtils: Curve
import Plots: mm

export apply_style!, logplot, draw!, save_fig

"""
    apply_style!(; title=16, guide=14, tick=12, legend=12)

Set the global font sizes. Call once at the top of `main`. The defaults are the
repo-wide values; pass a keyword only where a figure genuinely needs to differ.
"""
apply_style!(; title::Int = 16, guide::Int = 14, tick::Int = 12, legend::Int = 12) =
    default(titlefont = font(title), guidefont = font(guide),
            tickfont  = font(tick),  legendfont = font(legend))

"""
    logplot(; ylabel, decades=0:6, kwargs...) -> Plots.Plot

Empty panel with logarithmic time on x and ticks at `10^k`. Any keyword
overrides the default of the same name, so `logplot(ylabel=…, yscale=:log10)`
gives a log-log panel without a second function.
"""
function logplot(; decades = 0:6, kwargs...)
    base = (; xlabel        = L"t",
              xscale        = :log10,
              xticks        = (10.0 .^ decades,
                               [latexstring("10^{$i}") for i in decades]),
              legend        = :bottomright,
              framestyle    = :box,
              grid          = false,
              size          = (800, 600),
              left_margin   = 6mm, right_margin  = 4mm,
              top_margin    = 2mm, bottom_margin = 6mm)
    plot(; merge(base, values(kwargs))...)
end

"""
    draw!(p, curves) -> p

Add every `Curve` to `p`, in order. The only `plot!` call in the codebase.

Warns when two curves end up with the same colour *and* linestyle, which is what
a cyclic palette does once a sweep outgrows it — the curves are then
indistinguishable in the legend and the figure is quietly wrong.
"""
function draw!(p, curves)
    seen = Dict{Tuple{Any,Symbol},AbstractString}()
    for c in curves
        prev = get(seen, (c.color, c.linestyle), nothing)
        isnothing(prev) ? (seen[(c.color, c.linestyle)] = c.label) :
            @warn "Indistinguishable curves: same colour and linestyle" first=prev second=c.label
        plot!(p, c.x, c.y; label = c.label, color = c.color,
              linestyle = c.linestyle, lw = c.linewidth, alpha = 0.85)
    end
    p
end

"""
    save_fig(fig, dir, stem; exts=(".pdf", ".png"))

Write `dir/stem.ext` for each extension, creating `dir`. PNG travels alongside
the PDF so a figure can be eyeballed without opening a viewer.
"""
function save_fig(fig, dir::AbstractString, stem::AbstractString;
                  exts = (".pdf", ".png"))
    mkpath(dir)
    for e in exts
        path = joinpath(dir, stem * e)
        savefig(fig, path)
        println("Saved: $path")
    end
end

end # module PlotStyle
