module Config

using TOML, Dates, JLD2

export ExperimentConfig, PlotConfig,
       load_experiment_config, load_plot_config,
       config_outdir, find_latest_result,
       save_with_metadata

# ==========================================================================
# EXPERIMENT CONFIG — used by compute scripts
# ==========================================================================

"""
Holds all parameters that fully describe one simulation experiment.
Loaded from a TOML file; nothing is hardcoded in the Julia source.
"""
struct ExperimentConfig
    # ── Meta ──────────────────────────────────────────────────────────────
    name::String
    description::String
    # ── Physics ───────────────────────────────────────────────────────────
    N::Int
    boundary::Symbol            # :fixed | :periodic
    system_type::Symbol         # :springs | :masses
    nonlinear::Symbol           # :alpha | :beta
    param_values::Vector{Float64}
    delta_values::Vector{Float64}
    initial_condition::Symbol   # :low | :high | :mode
    initial_energy::Float64
    init_mode::Union{Int,Nothing}  # explicit mode index; nothing → derive from initial_condition
    # ── Simulation ────────────────────────────────────────────────────────
    TMAX::Float64
    T_block::Float64
    DT::Float64
    save_every::Int
    downsample::Int
    debug::Bool
    # ── Output ────────────────────────────────────────────────────────────
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

# ==========================================================================
# PLOT CONFIG — used by plotting scripts
# ==========================================================================

"""
Points a plotting script to a specific .jld2 result file
and provides all plot-level parameters.
"""
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
    color_scale::Symbol         # :log | :linear
    # ── Filters (empty = include all) ─────────────────────────────────────
    filter_params::Vector{Float64}
    filter_deltas::Vector{Float64}
end

function load_plot_config(path::String)::PlotConfig
    isfile(path) || error("Config file not found: $path")
    d = TOML.parsefile(path)

    data     = d["data"]
    analysis = get(d, "analysis", Dict{String,Any}())
    plot     = get(d, "plot",     Dict{String,Any}())
    out      = get(d, "output",   Dict{String,Any}())

    raw_t_max = Float64(get(plot, "t_max", 1e18))
    t_max_val = raw_t_max >= 1e17 ? Inf : raw_t_max
    PlotConfig(
        data["input_file"],
        Float64(get(analysis, "smooth_delta", 0.6)),
        get(out, "entropy_dir",  "results/figures/entropy"),
        get(out, "xi_dir",       "results/figures/xi"),
        get(out, "heatmaps_dir", "results/figures/heatmaps"),
        Float64(get(plot, "t_min", 0.0)),
        t_max_val,
        Int(get(plot, "time_subsample", 1)),
        Int(get(plot, "max_delta_cols", 3)),
        Symbol(get(plot, "color_scale", "log")),
        Float64.(get(get(d, "filter", Dict{String,Any}()), "param_values", Float64[])),
        Float64.(get(get(d, "filter", Dict{String,Any}()), "delta_values", Float64[])),
    )
end

# ==========================================================================
# UTILITIES
# ==========================================================================

"""Deterministic output directory for an experiment: base_dir/name/"""
config_outdir(cfg::ExperimentConfig) = joinpath(cfg.base_dir, cfg.name)

"""Returns the path of the most-recently-modified .jld2 inside a directory."""
function find_latest_result(experiment_dir::String)::String
    files = filter(f -> endswith(f, ".jld2"), readdir(experiment_dir; join=true))
    isempty(files) && error("No .jld2 files found in: $experiment_dir")
    sort(files; by=mtime)[end]
end

"""
Saves `results` + full reproducibility metadata to `outpath` (.jld2),
and copies the TOML config alongside it for human inspection.

Metadata fields written:
  config_toml, git_sha, julia_ver, project_toml, timestamp, host, pbs_jobid
"""
function save_with_metadata(outpath::String, config_path::String; kwargs...)
    git_sha    = try readchomp(`git rev-parse --short HEAD`) catch _ "" end
    julia_ver  = string(VERSION)
    project    = try read("Project.toml", String) catch _ "" end
    ts         = string(Dates.now())
    host       = get(ENV, "HOSTNAME", get(ENV, "HOST", "unknown"))
    jobid      = get(ENV, "PBS_JOBID", "")
    cfg_text   = read(config_path, String)

    jldsave(outpath;
        config_toml = cfg_text,
        git_sha     = git_sha,
        julia_ver   = julia_ver,
        project_toml = project,
        timestamp   = ts,
        host        = host,
        pbs_jobid   = jobid,
        kwargs...
    )

    # Copy config next to results for human inspection
    cp(config_path, joinpath(dirname(outpath), "config.toml"); force=true)
    println("Metadata saved: git=$(git_sha), host=$(host), t=$(ts)")
end

end
