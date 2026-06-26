module Config

using TOML

export ExperimentConfig, PlotConfig,
       load_experiment_config, load_plot_config,
       default_plot_config

struct ExperimentConfig
    name::String
    description::String
    N::Int
    boundary::Symbol
    system_type::Symbol
    nonlinear::Symbol
    param_values::Vector{Float64}
    delta_values::Vector{Float64}
    initial_condition::Symbol
    initial_energy::Float64
    init_mode::Union{Int,Nothing}
    TMAX::Float64
    T_block::Float64
    DT::Float64
    save_every::Int
    downsample::Int
    debug::Bool
    base_dir::String
end

function load_experiment_config(path::String)::ExperimentConfig
    isfile(path) || error("Config file not found: $path")
    d = TOML.parsefile(path)

    m = d["experiment"]
    p = d["physics"]
    s = d["simulation"]
    o = d["output"]

    ExperimentConfig(
        m["name"],
        get(m, "description", ""),
        Int(p["N"]),
        Symbol(p["boundary"]),
        Symbol(p["system_type"]),
        Symbol(p["nonlinear"]),
        Float64.(p["param_values"]),
        Float64.(p["delta_values"]),
        Symbol(get(p, "initial_condition", "low")),
        Float64(get(p, "initial_energy", 0.45)),
        haskey(p, "init_mode") ? Int(p["init_mode"]) : nothing,
        Float64(s["TMAX"]),
        Float64(s["T_block"]),
        Float64(get(s, "DT", 0.05)),
        Int(get(s, "save_every", 1000)),
        Int(get(s, "downsample", 2)),
        Bool(get(s, "debug", true)),
        o["base_dir"],
    )
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
        #:log,
        :linear,
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
