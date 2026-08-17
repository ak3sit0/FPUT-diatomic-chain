"""
    plot_heatmap_grid.jl

Render a grid of modal-energy heatmaps for parameter sweeps using the same
layout and style as the recovery script, but wired to this repository's data
and config format.

Usage:
  julia --project=. scripts/plot_heatmap_grid.jl <plot_config.toml | results.jld2>
"""

using CairoMakie, JLD2, LaTeXStrings, Colors
include("../../src/config.jl"); using Main.Config

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
    config = data["config"]
    if isa(config, String) && isfile(config)
        config = Config.load_experiment_config(config)
    end

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
    z_raw = Float64.(result.modal_E)
    if size(z_raw, 1) > size(z_raw, 2)
        z_raw = z_raw'
    end

    t_idx = findall(t -> cfg.t_min <= t <= cfg.t_max, t_raw)
    t_vals = t_raw[t_idx][1:cfg.time_subsample:end]
    z_vals = z_raw[:, t_idx][:, 1:cfg.time_subsample:end]
    z_plot = clamp.(z_vals, cfg.clamp_min, cfg.clamp_max)

    # If the plotting uses a log-scale x axis, replace zeros (log undefined).
    if any(t_vals .<= 0)
        pos = t_vals[t_vals .> 0]
        if !isempty(pos)
            minpos = minimum(pos)
            t_vals = [t <= 0 ? minpos/10 : t for t in t_vals]
        else
            # fallback tiny positive number if no positive times available
            t_vals = [t <= 0 ? 1e-12 : t for t in t_vals]
        end
    end

    PanelData(t_vals, 1:size(z_raw, 1), z_plot)
end

function find_result(dataset::HeatmapDataset, param, delta)
    idx = findfirst(r -> isapprox(r.param, param; rtol=1e-6) && isapprox(r.Delta, delta; rtol=1e-6), dataset.results)
    isnothing(idx) ? nothing : dataset.results[idx]
end

function make_panel_axis(ga, row, col, panel::PanelData, dataset::HeatmapDataset, cfg)
    n_rows = length(dataset.param_values)
    delta = dataset.delta_values[col]
    m = panel.mode_range
    y_tick_coords = range(first(m), last(m), length=5)
    y_tick_labels = [string(round(Int, y)) for y in y_tick_coords]

    # Create axis; if we have positive times use log10 scale with power-of-10 ticks
    pos = panel.t_vals[panel.t_vals .> 0]
    if !isempty(pos)
        minpos = minimum(pos)
        # prefer ticks at integer powers of ten starting at 10^0 up to max
        max_exp = floor(Int, log10(panel.t_vals[end]))
        candidate_exps = 0:max_exp
        xt_all = 10.0 .^ candidate_exps
        # keep only ticks inside the current panel time range
        xt = xt_all[(xt_all .>= panel.t_vals[1]) .& (xt_all .<= panel.t_vals[end])]
        xlabels = [string(Int(x)) for x in xt]
        # fallback: if no powers-of-ten fall inside range, revert to full computed range
        if isempty(xt)
            exp_min = floor(Int, log10(minpos))
            exp_max = ceil(Int, log10(panel.t_vals[end]))
            exps = exp_min:exp_max
            xt = 10.0 .^ exps
            xlabels = [string(Int(x)) for x in xt]
        end

        ax = Axis(ga[row, col];
            title              = row == 1 ? latexstring("\\Delta \\kappa = $delta") : "",
            titlesize          = cfg.titlesize,
            titlefont          = cfg.font,
            xticklabelsize     = cfg.ticksize,
            xticklabelfont     = cfg.font,
            yticklabelsize     = cfg.ticksize,
            yticklabelfont     = cfg.font,
            xticklabelsvisible = row == n_rows,
            yticklabelsvisible = col == 1,
            xscale             = log10,
            yticks             = (y_tick_coords, y_tick_labels),
            xticks             = (xt, xlabels),
        )
    else
        ax = Axis(ga[row, col];
            title              = row == 1 ? latexstring("\\Delta \\kappa = $delta") : "",
            titlesize          = cfg.titlesize,
            titlefont          = cfg.font,
            xticklabelsize     = cfg.ticksize,
            xticklabelfont     = cfg.font,
            yticklabelsize     = cfg.ticksize,
            yticklabelfont     = cfg.font,
            xticklabelsvisible = row == n_rows,
            yticklabelsvisible = col == 1,
            xscale             = identity,
            yticks             = (y_tick_coords, y_tick_labels),
        )
    end
    return ax
end

function render_heatmap_panel!(ax, panel::PanelData, cfg)
    hm = heatmap!(ax, panel.t_vals, panel.mode_range, panel.z_plot';
        colormap = cgrad([:white, "#B2D9FF", "#5999F2", "#3359CC", "#0D4CB3"]),
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

    if cfg.color_scale == :log
        Colorbar(fig[1, 2], hm_ref;
            label = "Log10(Energy)",
            labelsize = cfg.labelsize,
            labelfont = cfg.font,
            ticklabelsize = cfg.ticksize,
            ticklabelfont = cfg.font,
            ticks = [1e-6, 1e-4, 1e-2, 1.0],
            tickformat = _ -> ["10^-6", "10^-4", "10^-2", "1"],
            width = 36,
        )
    else
        Colorbar(fig[1, 2], hm_ref;
            label = "Energy",
            labelsize = cfg.labelsize,
            labelfont = cfg.font,
            ticklabelsize = cfg.ticksize,
            ticklabelfont = cfg.font,
            width = 36,
        )
    end

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
