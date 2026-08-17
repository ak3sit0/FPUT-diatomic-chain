module PlottingUtils

using JLD2, Statistics
using ..FPUTAnalysis

export cyc,
       PALETTE_DELTA, PALETTE_CATEGORICAL, LINESTYLES, LINEWIDTHS,
       BLUES_STOPS, BLUE_DARK, BLUE_MID, BLUE_LIGHT, GRAY_GUIDE,
       HALO_BLUE, HALO_TEAL, HALO_RED,
       PALETTE_COUPLING, RESONANCE_LINE,
       SPEC_FBC, SPEC_PBC, SPEC_PBC_COMPLETE,
       FTMLE_PALETTE, FTMLE_LINESTYLES,
       SMOOTH_DELTA,
       entropy_series, thermalization_time, prepare_ts, logdownsample,
       load_results_by_delta

# Backend-agnostic on purpose: palettes are plain hex strings / Symbols so both
# Plots.jl and CairoMakie can consume them without this module depending on either.

# ── Cyclic indexing ──────────────────────────────────────────────────────────

"""
    cyc(v, i)

Cyclic element access: `v[mod1(i, length(v))]`. Use for every palette/style
lookup so adding a curve past the end of a palette wraps instead of throwing.
"""
@inline cyc(v, i::Integer) = v[mod1(i, length(v))]

# ── Palettes ─────────────────────────────────────────────────────────────────

"""Sequential blues for Δκ sweeps (dispersion, resonance level curves)."""
const PALETTE_DELTA = [:darkblue, :steelblue, :cornflowerblue, :deepskyblue, :lightblue]

const LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]
const LINEWIDTHS = [2.4, 2.9, 3.3, 3.6]

"""
Categorical palette for sweeps whose curves are unordered labels rather than a
magnitude (e.g. system size N), where the sequential blues would read as a trend.
"""
const PALETTE_CATEGORICAL = [:blue, :red, :green, :orange, :purple, :brown, :magenta]

"""Colormap stops for modal-energy heatmaps. Build with `cgrad(BLUES_STOPS)`."""
const BLUES_STOPS = [:white, "#B2D9FF", "#5999F2", "#3359CC", "#0D4CB3"]

"""Colormap for |Γ| coupling heatmaps."""
const PALETTE_COUPLING = :Blues

"""Contrast color for resonance-manifold contours drawn over `PALETTE_COUPLING`."""
const RESONANCE_LINE = :orangered

const BLUE_DARK  = "#0D4CB3"   # lines, markers
const BLUE_MID   = "#5999F2"   # ±σ bands
const BLUE_LIGHT = "#B2D9FF"
const GRAY_GUIDE = :gray       # reference hlines

"""(marker, halo) colour pairs for σ-as-halo scatter plots."""
const HALO_BLUE = (:darkblue, :steelblue)
const HALO_TEAL = (:steelblue, :lightblue)
const HALO_RED  = (:darkred, :salmon)

# ── Entropy-figure specs: (Δκ, color, linestyle, linewidth) ──────────────────

const SPEC_FBC = [
    (0.05, "#1a3a6b", :solid,   4.0),
    (0.1,  "#2171b5", :dash,    2.5),
    (0.5,  "#d62728", :dashdot, 2.0),
    (0.7,  "#f4845f", :dot,     1.5),
]

const SPEC_PBC = [
    (0.05, "#1a3a6b", :solid,   4.0),
    (0.1,  "#2171b5", :dash,    2.5),
    (0.2,  "#6baed6", :dot,     2.5),
    (0.5,  "#d62728", :dashdot, 2.0),
    (0.7,  "#f4845f", :solid,   1.5),
]

const SPEC_PBC_COMPLETE = [
    (0.05, "#1a3a6b", :solid,   4.0),
    (0.1,  "#2171b5", :dash,    2.5),
    (0.2,  "#6baed6", :dot,     2.5),
    (0.3,  "#b3d9ff", :dashdot, 2.0),
    (0.5,  "#d62728", :dashdot, 2.0),
    (0.7,  "#f4845f", :solid,   1.5),
    (0.9,  "#8b0000", :dot,     1.5),
]

# ── ftMLE palette (contrasting blues → teal → violet) ────────────────────────

const FTMLE_PALETTE = ["#3399DB", "#2178B5", "#1AA6C2", "#5938A6", "#0D8080"]
const FTMLE_LINESTYLES = [:dot, :dash, :dashdot, :dashdotdot, :solid]

# ── Shared analysis helpers ──────────────────────────────────────────────────

"""Default growing-window smoothing factor Δ for S̄(t)."""
const SMOOTH_DELTA = 0.6

"""
    entropy_series(modal_E; delta=SMOOTH_DELTA) -> Vector{Float64}

Spectral entropy S̄(t) of an `N × nt` modal-energy matrix. Thin wrapper over
`FPUTAnalysis.spectral_entropy` so plotting scripts share the compute path with
the analysis module instead of reimplementing it.
"""
entropy_series(modal_E::AbstractMatrix; delta::Float64=SMOOTH_DELTA) =
    FPUTAnalysis.spectral_entropy(Matrix{Float64}(modal_E), delta)

"""
    thermalization_time(t, E_opt; threshold=0.9) -> Float64

First time at which `E_opt` reaches `threshold` of its asymptotic value (the mean
over the last 10% of the series). `NaN` when it never does.

Note this is the **optical-energy** definition used by the figure scripts. It is
*not* the same quantity as the `T_therm_*` fields stored by
`compute_ensemble.jl`, which threshold the normalized entropy instead — see
`docs/physics_diagnostics.md`.
"""
function thermalization_time(t::AbstractVector, E_opt::AbstractVector; threshold::Float64=0.9)
    (isempty(E_opt) || length(t) != length(E_opt)) && return NaN
    E_inf = mean(@view E_opt[max(1, round(Int, 0.9 * length(E_opt))):end])
    idx = findfirst(>=(threshold * E_inf), E_opt)
    isnothing(idx) ? NaN : t[idx]
end

"""
    prepare_ts(res) -> (t, E)

Extract `(scaled_t, modal_E)` from a result NamedTuple as dense `Float64`,
oriented `N × nt`, keeping only strictly positive times (log-scale x axis).
"""
function prepare_ts(res)
    t = Float64.(Vector(res.scaled_t))
    E = Float64.(Matrix(res.modal_E))
    size(E, 1) > size(E, 2) && (E = E')
    idx = findall(>(1e-3), t)
    t[idx], E[:, idx]
end

"""
    logdownsample(t, y, n_max=3000) -> (t, y)

Keep at most one sample per logarithmic bin, so log-x figures stay light
without visibly changing the curve.
"""
function logdownsample(t, y, n_max::Int=3000)
    length(t) <= n_max && return t, y
    log_edges = range(log10(t[1]), log10(t[end]); length=n_max + 1)
    idx = Int[]
    for i in 1:n_max
        lo, hi = 10^log_edges[i], 10^log_edges[i+1]
        j = findfirst(x -> lo <= x < hi, t)
        isnothing(j) || push!(idx, j)
    end
    isempty(idx) && return t, y
    t[idx], y[idx]
end

"""
    load_results_by_delta(paths...) -> Dict{Float64,Any}

Index results by Δκ across one or more JLD2 files. Earlier paths win on
collision, so a base sweep can be extended by later files without overriding it.
"""
function load_results_by_delta(paths::AbstractString...)
    by_delta = Dict{Float64, Any}()
    for path in paths
        for res in JLD2.load(path)["results"]
            d = Float64(res.Delta)
            haskey(by_delta, d) || (by_delta[d] = res)
        end
    end
    by_delta
end

end # module PlottingUtils
