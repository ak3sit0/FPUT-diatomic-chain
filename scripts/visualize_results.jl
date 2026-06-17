"""
    visualize_results.jl

Lightweight visualization for sweep results using CairoMakie for heatmap grids
and Plots for time-series. Intended to reproduce the multi-panel figures from
your recovery scripts for thesis-quality output.

Usage:
  julia --project=. scripts/visualize_results.jl results/raw/sweep_results_YYYY-MM-DD.jld2
"""

using JLD2, CairoMakie, LaTeXStrings
include("../src/fput_analysis.jl"); using .FPUTAnalysis
using Colors

# Match the recovery script palette and style
const PALETTE = cgrad([:white, "#B2D9FF", "#5999F2", "#3359CC", "#0D4CB3"])
const FONT = "Times New Roman"
const TITLESIZE = 35
const LABELSIZE = 45
const TICKSIZE = 30

function load_sweep(path::String)
    data = load(path)
    return data["results"], data["config"]
end

function plot_grid(results, outpath::String)
    # Recreate the recovery grid layout: rows=params, cols=deltas
    mkpath("results/figures")

    # Build unique sorted params and deltas
    params = sort(unique(Float64[r.param for r in results]))
    deltas = sort(unique(Float64[r.Delta for r in results]))
    n_rows = length(params)
    n_cols = length(deltas)

    # Load configuration constants from this file
    cfg = (
        colormap = PALETTE,
        clamp_min = 1e-6,
        clamp_max = 1.0,
        px_per_unit = 2,
        font = FONT,
        titlesize = TITLESIZE,
        labelsize = LABELSIZE,
        ticksize = TICKSIZE,
        cell_w = 600,
        cell_h = 500,
        margin_w = 200,
        margin_h = 150,
        gap = 38,
        color_scale = :log,
        t_min = 1.0,
    )

    fig = Figure(size = (cfg.cell_w * n_cols + cfg.margin_w, cfg.cell_h * n_rows + cfg.margin_h), font = cfg.font)
    ga = fig[1,1] = GridLayout()

    hm_ref = nothing

    for i in 1:n_rows
        # rows go top-to-bottom; flip params so largest is on top
        param = params[n_rows - i + 1]
        for j in 1:n_cols
            delta = deltas[j]
            # find result matching param & delta
            idx = findfirst(r -> isapprox(Float64(r.param), param, rtol=1e-8) && isapprox(Float64(r.Delta), delta, rtol=1e-8), results)
            if isnothing(idx)
                continue
            end
            r = results[idx]
            # Prepare panel data
            t_raw = Vector{Float64}(r.scaled_t)
            z_raw = Matrix{Float64}(r.modal_E)
            if size(z_raw,1) > size(z_raw,2)
                z_raw = z_raw'
            end
            # Time filtering to avoid log10(0) on axis
            t_idx = findall(t -> t >= cfg.t_min, t_raw)
            if isempty(t_idx)
                t_vals = t_raw
                z_subset = z_raw
            else
                t_vals = t_raw[t_idx]
                z_subset = z_raw[:, t_idx]
            end
            z_plot = clamp.(z_subset, cfg.clamp_min, cfg.clamp_max)

            ax = Axis(ga[i, j];
                titlesize = cfg.titlesize,
                xticklabelsize = cfg.ticksize,
                yticklabelsize = cfg.ticksize,
                xscale = log10,
            )

            hm = heatmap!(ax, t_vals, 1:size(z_plot,1), log10.(z_plot .+ 1e-15)'; colormap = cfg.colormap)
            isnothing(hm_ref) && (hm_ref = hm)
        end
    end

    # Add labels and colorbar (LaTeX labels)
    Label(fig[1, 1, Bottom()], L"t"; fontsize = cfg.labelsize, padding = (0,0,0,40), font = cfg.font)
    Label(fig[1, 1, Left()], L"Mode\ index"; fontsize = cfg.labelsize, rotation = pi/2, padding = (0,60,0,0), font = cfg.font)
    cbticks_vals = Float64[-6.0, -4.0, -2.0, 0.0]
    cbticks_lbls = [L"10^{-6}", L"10^{-4}", L"10^{-2}", L"1"]
    Colorbar(fig[1, 2], hm_ref; label = L"Energy", labelsize = cfg.labelsize, ticklabelsize=cfg.ticksize, width=36, ticks=(cbticks_vals, cbticks_lbls))

    save(joinpath("results/figures", "heatmap_grid_fromsweep.pdf"), fig)
    println("Saved results/figures/heatmap_grid_fromsweep.pdf")
    println("Saved results/figures/heatmap_grid_fromsweep.pdf")
end

function main()
    if isempty(ARGS)
        println("Usage: julia scripts/visualize_results.jl <sweep_jld2>")
        return
    end
    results, cfg = load_sweep(ARGS[1])
    plot_grid(results, ARGS[1])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
