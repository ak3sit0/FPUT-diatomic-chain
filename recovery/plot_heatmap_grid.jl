# Load precomputed heatmap data and render a (params × deltas) grid of energy heatmaps.
# Usage: julia --project=. examples/plot_heatmap_grid.jl configs/plot_heatmap_alpha_periodic.toml

using CairoMakie, JLD2, FileIO, LaTeXStrings
include("../src/config.jl"); using Main.Config

# ==========================================================================
# CONFIGURATION
# ==========================================================================

function build_plot_config(config_path::Union{String,Nothing} = nothing)
    # Load TOML values when a config file is provided,
    # then merge with (or fall back to) hard-coded defaults.
    pcfg = isnothing(config_path) ? nothing : Config.load_plot_config(config_path)

    default_file = "results/data/modal_energies/springs_alpha_periodic_delta_k_00103/" *
                   "modal_energies_springs_alpha_periodic_2026-02-06.jld2"
    (
        datafile       = isnothing(pcfg) ? default_file : pcfg.input_file,
        outdir         = isnothing(pcfg) ? "results/figures/heatmaps" : pcfg.outdir_heatmaps,
        t_min          = isnothing(pcfg) ? 1.0   : pcfg.t_min,
        t_max          = isnothing(pcfg) ? Inf   : pcfg.t_max,
        time_subsample = isnothing(pcfg) ? 1     : pcfg.time_subsample,
        max_delta_cols = isnothing(pcfg) ? 3     : pcfg.max_delta_cols,
        color_scale    = isnothing(pcfg) ? :log  : pcfg.color_scale,
        filter_params  = isnothing(pcfg) ? Float64[] : pcfg.filter_params,
        filter_deltas  = isnothing(pcfg) ? Float64[] : pcfg.filter_deltas,
        clamp_min      = 1e-6,
        clamp_max      = 1.0,
        colormap       = cgrad([:white, "#B2D9FF", "#5999F2", "#3359CC", "#0D4CB3"]),
        px_per_unit    = 2,           # 1 ≈ 72 DPI, 2 ≈ 144 DPI
        font           = "Times New Roman",
        titlesize      = 35,
        labelsize      = 45,
        ticksize       = 30,
        cell_w         = 600,
        cell_h         = 500,
        margin_w       = 200,
        margin_h       = 150,
        gap            = 38,
    )
end

# ==========================================================================
# DATA LOADING
# ==========================================================================

struct HeatmapDataset
    results::Vector
    config::Any
    param_values::Vector{Float64}   # rows (nonlinearity), sorted ascending
    delta_values::Vector{Float64}   # columns (Δκ or Δm), sorted ascending
end

function load_heatmap_data(datafile::String, max_cols::Real = Inf;
                            filter_params::Vector{Float64}=Float64[],
                            filter_deltas::Vector{Float64}=Float64[])
    isfile(datafile) || error("File not found: $datafile")
    println("Loading: $datafile")
    data    = FileIO.load(datafile)
    results = data["results"]
    config  = data["config"]

    # Apply filters (empty = no filter)
    if !isempty(filter_params)
        results = filter(r -> any(p -> isapprox(Float64(r.param), p, rtol=1e-6), filter_params), results)
    end
    if !isempty(filter_deltas)
        results = filter(r -> any(d -> isapprox(Float64(r.Delta), d, rtol=1e-6), filter_deltas), results)
    end

    params = sort(unique(Float64[r.param for r in results]))
    deltas = sort(unique(Float64[r.Delta for r in results]))

    n_cols = min(length(deltas), Int(floor(max_cols)))
    deltas = deltas[1:n_cols]

    HeatmapDataset(results, config, params, deltas)
end

# ==========================================================================
# PANEL DATA PREPARATION  (pure — no plotting)
# ==========================================================================

struct PanelData
    t_vals::Vector{Float64}
    mode_range::UnitRange{Int}
    z_plot::Matrix{Float64}    # (modes × time_steps), ready to render
end

function prepare_panel_data(result, cfg)
    t_raw = Vector{Float64}(result.scaled_t)
    z_raw = Matrix{Float64}(result.modal_E)

    # Guarantee (modes × time)
    if size(z_raw, 1) > size(z_raw, 2)
        z_raw = z_raw'
    end

    # Time range filter
    t_idx  = findall(t -> cfg.t_min <= t <= cfg.t_max, t_raw)
    t_vals = t_raw[t_idx][1:cfg.time_subsample:end]
    z_vals = z_raw[:, t_idx][:, 1:cfg.time_subsample:end]

    # Clamp for color scale stability
    z_plot = clamp.(z_vals, cfg.clamp_min, cfg.clamp_max)

    PanelData(t_vals, 1:size(z_raw, 1), z_plot)
end

function find_result(dataset::HeatmapDataset, param, delta)
    idx = findfirst(r -> r.param ≈ param && r.Delta ≈ delta, dataset.results)
    isnothing(idx) ? nothing : dataset.results[idx]
end

# ==========================================================================
# AXIS & PANEL RENDERING
# ==========================================================================

function make_panel_axis(ga, row, col, panel::PanelData, dataset::HeatmapDataset, cfg)
    n_rows = length(dataset.param_values)
    delta  = dataset.delta_values[col]
    m      = panel.mode_range
    y_tick_coords  = range(first(m), last(m), length=5)
    y_tick_labels  = [string(round(Int, y)) for y in y_tick_coords]

    Axis(ga[row, col];
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
    )
end

function render_heatmap_panel!(ax, panel::PanelData, cfg)
    scale  = cfg.color_scale == :log ? log10 : identity
    crange = cfg.color_scale == :log ? (cfg.clamp_min, cfg.clamp_max) : (0.0, cfg.clamp_max)

    hm = heatmap!(ax, panel.t_vals, panel.mode_range, panel.z_plot';
        colormap   = cfg.colormap,
        colorscale = scale,
        colorrange = crange,
    )
    xlims!(ax, panel.t_vals[1], panel.t_vals[end])
    ylims!(ax, first(panel.mode_range) - 0.5, last(panel.mode_range) + 0.5)
    hm
end

# ==========================================================================
# FIGURE ASSEMBLY
# ==========================================================================

function add_labels_and_colorbar!(fig, ga, hm_ref, cfg)
    Label(fig[1, 1, Bottom()], "Scaled Time";
        fontsize = cfg.labelsize, padding = (0, 0, 0, 40), font = cfg.font)
    Label(fig[1, 1, Left()], "Mode Index";
        fontsize = cfg.labelsize, rotation = pi/2, padding = (0, 60, 0, 0), font = cfg.font)

    if cfg.color_scale == :log
        Colorbar(fig[1, 2], hm_ref;
            label         = "Energy (log10)",
            labelsize     = cfg.labelsize,
            labelfont     = cfg.font,
            ticklabelsize = cfg.ticksize,
            ticklabelfont = cfg.font,
            ticks         = [1e-6, 1e-4, 1e-2, 1.0],
            tickformat    = _ -> ["10^-6", "10^-4", "10^-2", "1"],
            width         = 36,        # make the bar thicker for readability
        )
    else
        Colorbar(fig[1, 2], hm_ref;
            label         = "Energy",
            labelsize     = cfg.labelsize,
            labelfont     = cfg.font,
            ticklabelsize = cfg.ticksize,
            ticklabelfont = cfg.font,
            width         = 36,
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
    ga  = fig[1, 1] = GridLayout()

    hm_ref = nothing

    for i in 1:n_rows
        # Rows go top-to-bottom; flip so the largest param value sits on top
        param = dataset.param_values[n_rows - i + 1]

        for j in 1:n_cols
            delta  = dataset.delta_values[j]
            result = find_result(dataset, param, delta)
            isnothing(result) && continue

            panel  = prepare_panel_data(result, cfg)
            ax     = make_panel_axis(ga, i, j, panel, dataset, cfg)
            hm     = render_heatmap_panel!(ax, panel, cfg)
            isnothing(hm_ref) && (hm_ref = hm)
        end
    end

    add_labels_and_colorbar!(fig, ga, hm_ref, cfg)
    return fig
end

# ==========================================================================
# SAVING
# ==========================================================================

function save_figure(fig, dataset::HeatmapDataset, cfg)
    mkpath(cfg.outdir)
    nonlinear = String(dataset.config.nonlinear)
    boundary  = String(dataset.config.boundary)
    fname     = joinpath(cfg.outdir, "heatmap_grid_$(nonlinear)_$(boundary)_makie.png")
    save(fname, fig; px_per_unit = cfg.px_per_unit)
    println("Saved: $fname")
    fname
end

# ==========================================================================
# MAIN ENTRY
# ==========================================================================

function main()
    config_path = isempty(ARGS) ? nothing : ARGS[1]
    if isnothing(config_path)
        println("Tip: pass a TOML config file for full reproducibility.\n" *
                "Usage: julia --project=. examples/plot_heatmap_grid.jl configs/plot_heatmap_alpha_periodic.toml")
    end
    cfg     = build_plot_config(config_path)
    dataset = load_heatmap_data(cfg.datafile, cfg.max_delta_cols;
                                filter_params=cfg.filter_params,
                                filter_deltas=cfg.filter_deltas)
    fig     = build_heatmap_figure(dataset, cfg)
    save_figure(fig, dataset, cfg)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
