module Config

using TOML

export ExperimentConfig, PlotConfig,
       load_experiment_config, load_plot_config,
       default_plot_config, resolve_config_path, as_experiment_config

struct ExperimentConfig
    system_type::Symbol
    boundary::Symbol
    nonlinear::Symbol
end

function load_experiment_config(path::String)::ExperimentConfig
    isfile(path) || error("Config file not found: $path")
    p = TOML.parsefile(path)["physics"]

    ExperimentConfig(
        Symbol(p["system_type"]),
        Symbol(p["boundary"]),
        Symbol(p["nonlinear"]),
    )
end

const CONFIG_DIRS = ("configs/production", "configs/hpc_templates", "configs/tests", "configs")

"""
    resolve_config_path(path) -> String | nothing

Locate a config TOML whose recorded path may be stale. Results saved before the
`configs/cases/` → `configs/production/` reorganization store paths like
`configs/cases/periodic_N64_production.toml`; this falls back to a basename
lookup (with and without the dropped `_production` suffix) across `configs/`.
Returns `nothing` when no candidate exists.
"""
function resolve_config_path(path::AbstractString)
    isfile(path) && return String(path)
    base = basename(path)
    candidates = unique([base, replace(base, "_production.toml" => ".toml")])
    for dir in CONFIG_DIRS, c in candidates
        p = joinpath(dir, c)
        isfile(p) && return p
    end
    nothing
end

"""
    as_experiment_config(config; fallback=ExperimentConfig(:springs, :periodic, :alpha))

Normalize the `config` entry stored in a results JLD2 — either an
`ExperimentConfig` or a (possibly stale) path string — into an
`ExperimentConfig`. Warns and returns `fallback` when the path cannot be
resolved, so plotting degrades to a default label instead of erroring on a
`String` that has no fields.
"""
function as_experiment_config(config;
                              fallback::ExperimentConfig = ExperimentConfig(:springs, :periodic, :alpha))
    config isa ExperimentConfig && return config
    if config isa AbstractString
        resolved = resolve_config_path(config)
        isnothing(resolved) || return load_experiment_config(resolved)
        @warn "Config referenced by the results file was not found; using defaults for labels/filenames" path=config fallback
    end
    fallback
end

struct PlotConfig
    input_file::String
    smooth_delta::Float64
    outdir_entropy::String
    outdir_xi::String
    outdir_heatmaps::String
    t_min::Float64
    t_max::Float64
    time_subsample::Int
    max_delta_cols::Int
    color_scale::Symbol
    filter_params::Vector{Float64}
    filter_deltas::Vector{Float64}
    clamp_min::Float64
    clamp_max::Float64
    px_per_unit::Int
    font::String
    titlesize::Int
    labelsize::Int
    ticksize::Int
    cell_w::Int
    cell_h::Int
    margin_w::Int
    margin_h::Int
    gap::Int
end

function default_plot_config(input_file::String)::PlotConfig
    PlotConfig(
        input_file,
        0.6,
        "results/figures/entropy",
        "results/figures/xi",
        "results/figures/heatmaps",
        1.0,
        Inf,
        1,
        3,
        :log,
        #:linear,
        Float64[],
        Float64[],
        1e-6,
        1.0,
        2,
        "Times New Roman",
        35,
        45,
        30,
        600,
        500,
        200,
        150,
        38,
    )
end

function load_plot_config(path::String)::PlotConfig
    isfile(path) || error("Config file not found: $path")
    d = TOML.parsefile(path)

    data = d["data"]
    analysis = get(d, "analysis", Dict{String, Any}())
    plot = get(d, "plot", Dict{String, Any}())
    out = get(d, "output", Dict{String, Any}())
    filter = get(d, "filter", Dict{String, Any}())
    heatmap = get(d, "heatmap", Dict{String, Any}())

    raw_t_max = Float64(get(plot, "t_max", 1e18))
    t_max_val = raw_t_max >= 1e17 ? Inf : raw_t_max

    PlotConfig(
        data["input_file"],
        Float64(get(analysis, "smooth_delta", 0.6)),
        get(out, "entropy_dir", "results/figures/entropy"),
        get(out, "xi_dir", "results/figures/xi"),
        get(out, "heatmaps_dir", get(out, "outdir_heatmaps", "results/figures/heatmaps")),
        Float64(get(plot, "t_min", 0.0)),
        t_max_val,
        Int(get(plot, "time_subsample", 1)),
        Int(get(plot, "max_delta_cols", 3)),
        Symbol(get(plot, "color_scale", get(heatmap, "scale", "log"))),
        Float64.(get(filter, "param_values", Float64[])),
        Float64.(get(filter, "delta_values", Float64[])),
        Float64(get(heatmap, "clamp_min", 1e-6)),
        Float64(get(heatmap, "clamp_max", 1.0)),
        Int(get(heatmap, "px_per_unit", 2)),
        String(get(heatmap, "font", "Times New Roman")),
        Int(get(heatmap, "titlesize", 35)),
        Int(get(heatmap, "labelsize", 45)),
        Int(get(heatmap, "ticksize", 30)),
        Int(get(heatmap, "cell_w", 600)),
        Int(get(heatmap, "cell_h", 500)),
        Int(get(heatmap, "margin_w", 200)),
        Int(get(heatmap, "margin_h", 150)),
        Int(get(heatmap, "gap", 38)),
    )
end

end # module
