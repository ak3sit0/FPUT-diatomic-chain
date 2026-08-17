"""
    plot_heatmap_grid.jl

Render a grid of modal-energy heatmaps for parameter sweeps using the same
layout and style as the recovery script, but wired to this repository's data
and config format.

Usage:
  julia --project=. scripts/plot_heatmap_grid.jl <plot_config.toml | results.jld2>
"""

using CairoMakie, JLD2, LaTeXStrings, Colors, Statistics
include("../../src/config.jl");         using .Config
include("../../src/fput_analysis.jl");  using .FPUTAnalysis
include("../../src/plotting_utils.jl"); using .PlottingUtils

function build_plot_config(input_arg::Union{String,Nothing} = nothing)
    if isnothing(input_arg)
        error("Usage: julia scripts/plot_heatmap_grid.jl <plot_config.toml | results.jld2>")
    elseif endswith(lowercase(input_arg), ".jld2")
        return Config.default_plot_config(input_arg)
    else
        return Config.load_plot_config(input_arg)
    end
end

struct HeatmapDataset
    results::Vector
    config::Any
    param_values::Vector{Float64}
    delta_values::Vector{Float64}
end

function load_heatmap_data(datafile::String; filter_params::Vector{Float64}=Float64[], filter_deltas::Vector{Float64}=Float64[])
    isfile(datafile) || error("File not found: $datafile")
    println("Loading: $datafile")
    data = load(datafile)
    results = data["results"]
    config = Config.as_experiment_config(data["config"])

    if !isempty(filter_params)
        results = filter(r -> any(p -> isapprox(Float64(r.param), p, rtol=1e-6), filter_params), results)
    end
    if !isempty(filter_deltas)
        results = filter(r -> any(d -> isapprox(Float64(r.Delta), d, rtol=1e-6), filter_deltas), results)
    end

    params = sort(unique(Float64[r.param for r in results]))
    deltas = sort(unique(Float64[r.Delta for r in results]))
    HeatmapDataset(results, config, params, deltas)
end

struct PanelData
    t_vals::Vector{Float64}
    mode_range::UnitRange{Int}
    z_plot::Matrix{Float64}
end

function prepare_panel_data(result, cfg)
    t_raw = Float64.(result.scaled_t)
    z_raw = result.modal_E
    size(z_raw, 1) > size(z_raw, 2) && (z_raw = permutedims(z_raw))

    # Window to the time range, then subsample. Slice once to avoid multiple copies.
    idx = findall(t -> cfg.t_min <= t <= cfg.t_max, t_raw)
    cfg.time_subsample > 1 && (idx = idx[1:cfg.time_subsample:end])

    t_vals = t_raw[idx]
    z_plot = clamp.(Float64.(z_raw[:, idx]), cfg.clamp_min, cfg.clamp_max)

    # Log-scale x axis: push non-positive times a decade below the smallest positive one.
    if any(<=(0), t_vals)
        pos = filter(>(0), t_vals)
        floor_t = isempty(pos) ? 1e-12 : minimum(pos) / 10
        t_vals = [t <= 0 ? floor_t : t for t in t_vals]
    end

    PanelData(t_vals, 1:size(z_raw, 1), z_plot)
end

function find_result(dataset::HeatmapDataset, param, delta)
    idx = findfirst(r -> isapprox(r.param, param; rtol=1e-6) && isapprox(r.Delta, delta; rtol=1e-6), dataset.results)
    isnothing(idx) ? nothing : dataset.results[idx]
end

"""
    decade_ticks(t_vals) -> (positions, labels) | nothing

Ticks at integer powers of ten inside the panel's time range; widens to the
enclosing decades when none fall inside. `nothing` when there is no positive
time to place on a log axis.
"""
function decade_ticks(t_vals)
    (isempty(t_vals) || !any(>(0), t_vals)) && return nothing
    xt = [10.0^e for e in 0:floor(Int, log10(t_vals[end])) if t_vals[1] <= 10.0^e <= t_vals[end]]
    if isempty(xt)
        lo = floor(Int, log10(minimum(filter(>(0), t_vals))))
        hi = ceil(Int, log10(t_vals[end]))
        xt = [10.0^e for e in lo:hi]
    end
    xt, [x >= 1 ? string(Int(x)) : string(x) for x in xt]
end

function make_panel_axis(ga, row, col, panel::PanelData, dataset::HeatmapDataset, cfg)
    n_rows = length(dataset.param_values)
    m = panel.mode_range
    y_tick_coords = range(first(m), last(m), length=5)
    ticks = decade_ticks(panel.t_vals)

    Axis(ga[row, col];
        title              = row == 1 ? latexstring("\\Delta \\kappa = $(dataset.delta_values[col])") : "",
        titlesize          = cfg.titlesize,
        titlefont          = cfg.font,
        xticklabelsize     = cfg.ticksize,
        xticklabelfont     = cfg.font,
        yticklabelsize     = cfg.ticksize,
        yticklabelfont     = cfg.font,
        xticklabelsvisible = row == n_rows,
        yticklabelsvisible = col == 1,
        xscale             = isnothing(ticks) ? identity : log10,
        yticks             = (y_tick_coords, [string(round(Int, y)) for y in y_tick_coords]),
        # Makie treats `automatic` as "pick your own"; only override when we have decades.
        (isnothing(ticks) ? () : (xticks = ticks,))...,
    )
end

function render_heatmap_panel!(ax, panel::PanelData, cfg)
    hm = heatmap!(ax, panel.t_vals, panel.mode_range, panel.z_plot';
        colormap = cgrad(BLUES_STOPS),
        colorscale = cfg.color_scale == :log ? log10 : identity,
        colorrange = cfg.color_scale == :log ? (cfg.clamp_min, cfg.clamp_max) : (0.0, cfg.clamp_max),
        rasterize = 4,
    )
    xlims!(ax, panel.t_vals[1], panel.t_vals[end])
    ylims!(ax, first(panel.mode_range) - 0.5, last(panel.mode_range) + 0.5)
    hm
end

function add_labels_and_colorbar!(fig, ga, hm_ref, cfg)
    Label(fig[1, 1, Bottom()], "Scaled Time";
        fontsize = cfg.labelsize, padding = (0, 0, 0, 40), font = cfg.font)
    Label(fig[1, 1, Left()], "Mode Index";
        fontsize = cfg.labelsize, rotation = pi/2, padding = (0, 60, 0, 0), font = cfg.font)

    is_log = cfg.color_scale == :log
    Colorbar(fig[1, 2], hm_ref;
        label         = is_log ? "Log10(Energy)" : "Energy",
        labelsize     = cfg.labelsize,
        labelfont     = cfg.font,
        ticklabelsize = cfg.ticksize,
        ticklabelfont = cfg.font,
        width         = 36,
        (is_log ? (ticks = [1e-6, 1e-4, 1e-2, 1.0],
                   tickformat = _ -> ["10^-6", "10^-4", "10^-2", "1"]) : ())...,
    )

    colgap!(ga, cfg.gap)
    rowgap!(ga, cfg.gap)
end

function build_heatmap_figure(dataset::HeatmapDataset, cfg)
    n_rows = length(dataset.param_values)
    n_cols = length(dataset.delta_values)

    fig = Figure(size = (cfg.cell_w * n_cols + cfg.margin_w, cfg.cell_h * n_rows + cfg.margin_h),
                 font = cfg.font)
    ga = fig[1, 1] = GridLayout()

    hm_ref = nothing
    for i in 1:n_rows
        param = dataset.param_values[n_rows - i + 1]
        for j in 1:n_cols
            delta = dataset.delta_values[j]
            result = find_result(dataset, param, delta)
            isnothing(result) && continue

            panel = prepare_panel_data(result, cfg)
            ax = make_panel_axis(ga, i, j, panel, dataset, cfg)
            hm = render_heatmap_panel!(ax, panel, cfg)
            isnothing(hm_ref) && (hm_ref = hm)
        end
    end

    add_labels_and_colorbar!(fig, ga, hm_ref, cfg)
    return fig
end

function save_figure(fig, dataset::HeatmapDataset, cfg)
    mkpath(cfg.outdir_heatmaps)
    nonlinear = String(dataset.config.nonlinear)
    boundary = String(dataset.config.boundary)
    fname = joinpath(cfg.outdir_heatmaps, "heatmap_grid_$(nonlinear)_$(boundary)_makie.png")
    save(fname, fig; px_per_unit = cfg.px_per_unit)
    println("Saved: $fname")
    fname
end

function main()
    cfg = build_plot_config(isempty(ARGS) ? nothing : ARGS[1])
    dataset = load_heatmap_data(cfg.input_file;
                                filter_params = cfg.filter_params,
                                filter_deltas = cfg.filter_deltas)
    fig = build_heatmap_figure(dataset, cfg)
    save_figure(fig, dataset, cfg)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
