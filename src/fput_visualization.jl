module FPUTVisualization

using CairoMakie, LaTeXStrings

export build_heatmap_figure, save_heatmap_figure

const DEFAULT = (
    clamp_min = 1e-6,
    clamp_max = 1.0,
    colormap = cgrad([:white, "#B2D9FF", "#5999F2", "#3359CC", "#0D4CB3"]),
    font = "Times New Roman",
    titlesize = 35,
    labelsize = 45,
    ticksize = 30,
    cell_w = 600,
    cell_h = 500,
    margin_w = 200,
    margin_h = 150,
    gap = 38,
    color_scale = :log,
    t_min = 1.0,
    px_per_unit = 2,
)

struct PanelData
    t_vals::Vector{Float64}
    mode_range::UnitRange{Int}
    z_plot::Matrix{Float64}
end

function prepare_panel_data(result; cfg=DEFAULT)
    t_raw = Vector{Float64}(result.scaled_t)
    z_raw = Matrix{Float64}(result.modal_E)
    if size(z_raw,1) > size(z_raw,2)
        z_raw = z_raw'
    end
    t_idx = findall(t -> cfg.t_min <= t <= typemax(Float64), t_raw)
    t_vals = t_raw[t_idx]
    z_vals = z_raw[:, t_idx]
    z_plot = clamp.(z_vals, cfg.clamp_min, cfg.clamp_max)
    PanelData(t_vals, 1:size(z_raw,1), z_plot)
end

function make_panel_axis(ga, row, col, panel::PanelData, dataset_params, cfg)
    n_rows = length(dataset_params)
    delta = dataset_params[col]
    m = panel.mode_range
    y_tick_coords = range(first(m), last(m), length=5)
    y_tick_labels = [string(round(Int, y)) for y in y_tick_coords]

    Axis(ga[row, col];
        title = row == 1 ? latexstring("\\Delta \\kappa = $delta") : "",
        titlesize = cfg.titlesize,
        titlefont = cfg.font,
        xticklabelsize = cfg.ticksize,
        xticklabelfont = cfg.font,
        yticklabelsize = cfg.ticksize,
        yticklabelfont = cfg.font,
        xticklabelsvisible = row == n_rows,
        yticklabelsvisible = col == 1,
        xscale = log10,
        yticks = (y_tick_coords, y_tick_labels),
    )
end

function render_heatmap_panel!(ax, panel::PanelData, cfg)
    scale = cfg.color_scale == :log ? log10 : identity
    crange = cfg.color_scale == :log ? (cfg.clamp_min, cfg.clamp_max) : (0.0, cfg.clamp_max)
    hm = heatmap!(ax, panel.t_vals, panel.mode_range, panel.z_plot';
        colormap = cfg.colormap,
        colorscale = scale,
        colorrange = crange,
    )
    xlims!(ax, panel.t_vals[1], panel.t_vals[end])
    ylims!(ax, first(panel.mode_range) - 0.5, last(panel.mode_range) + 0.5)
    hm
end

function add_labels_and_colorbar!(fig, ga, hm_ref, cfg)
    Label(fig[1, 1, Bottom()], "Scaled Time"; fontsize = cfg.labelsize, padding = (0,0,0,40), font = cfg.font)
    Label(fig[1, 1, Left()], "Mode Index"; fontsize = cfg.labelsize, rotation = pi/2, padding = (0,60,0,0), font = cfg.font)

    if cfg.color_scale == :log
        Colorbar(fig[1, 2], hm_ref;
            label = "Energy (log10)",
            labelsize = cfg.labelsize,
            labelfont = cfg.font,
            ticklabelsize = cfg.ticksize,
            ticklabelfont = cfg.font,
            ticks = [1e-6, 1e-4, 1e-2, 1.0],
            width = 36,
        )
    else
        Colorbar(fig[1, 2], hm_ref; label = "Energy", labelsize = cfg.labelsize, labelfont = cfg.font, ticklabelsize = cfg.ticksize, width = 36)
    end

    # layout spacing
    rowgap!(ga, cfg.gap)
    colgap!(ga, cfg.gap)
end

function build_heatmap_figure(results; cfg=DEFAULT)
    params = sort(unique(Float64[r.param for r in results]))
    deltas = sort(unique(Float64[r.Delta for r in results]))
    n_rows = length(params)
    n_cols = length(deltas)

    fig = Figure(size = (cfg.cell_w * n_cols + cfg.margin_w, cfg.cell_h * n_rows + cfg.margin_h), font = cfg.font)
    ga = fig[1, 1] = GridLayout()
    hm_ref = nothing

    for i in 1:n_rows
        param = params[n_rows - i + 1]
        for j in 1:n_cols
            delta = deltas[j]
            idx = findfirst(r -> isapprox(Float64(r.param), param, rtol=1e-8) && isapprox(Float64(r.Delta), delta, rtol=1e-8), results)
            isnothing(idx) && continue
            r = results[idx]
            panel = prepare_panel_data(r; cfg=cfg)
            ax = make_panel_axis(ga, i, j, panel, deltas, cfg)
            hm = render_heatmap_panel!(ax, panel, cfg)
            isnothing(hm_ref) && (hm_ref = hm)
        end
    end

    add_labels_and_colorbar!(fig, ga, hm_ref, cfg)
    return fig
end

function save_heatmap_figure(fig, outdir::String; cfg=DEFAULT)
    mkpath(outdir)
    fname = joinpath(outdir, "heatmap_grid_refactored.png")
    save(fname, fig; px_per_unit = cfg.px_per_unit)
    return fname
end

end # module
