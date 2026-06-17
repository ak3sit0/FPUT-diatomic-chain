"""
examples/plot_entropy_N_timeseries.jl

Lee un archivo JLD2 generado por `compute_trajectories_Nsweep.jl` y plotea
múltiples curvas `S(t)` superpuestas, una por cada `N`.

Usage:
  julia --project=. examples/plot_entropy_N_timeseries.jl results/data/sweep_N_results_YYYY-MM-DD.jld2
"""

using JLD2, Plots, LaTeXStrings
gr()
include("../src/fput_analysis.jl"); using .FPUTAnalysis

const PLOT_MAX_POINTS = 2000 # max points per curve to plot (downsample if longer)

function main()
    if isempty(ARGS)
        println("Usage: julia examples/plot_entropy_N_timeseries.jl <results.jld2>")
        return
    end
    data = load(ARGS[1])
    results = data["results"]

    # helper to robustly extract fields from NamedTuple or Dict-like entries
    getval(r, name) = try
        getproperty(r, Symbol(name))
    catch
        try
            r[name]
        catch
            try
                r[string(name)]
            catch
                nothing
            end
        end
    end

    palette = [:blue, :red, :green, :orange, :purple, :brown, :magenta]
    p = plot(xlabel=L"cycles", ylabel=L"S(t)", title="Spectral Entropy vs time for different N",
             lw=2.5, legend=:outertopright, grid=true, size=(1000,600))

    for (i, r) in enumerate(results)
        t_raw = getval(r, "scaled_t")
        modal_E_raw = getval(r, "modal_E")
        S_raw = getval(r, "entropy")

        if t_raw === nothing
            println("Skipping entry $i: no time vector found")
            continue
        end

        t = Float64.(t_raw)

        # If entropy missing but modal_E present, compute it here
        if S_raw === nothing && modal_E_raw !== nothing
            modal_E = Array{Float64}(modal_E_raw)
            S = FPUTAnalysis.spectral_entropy(modal_E, 0.6)
        elseif S_raw !== nothing
            S = Float64.(S_raw)
        else
            println("Skipping entry $i: neither entropy nor modal_E available")
            continue
        end

        # Downsample for plotting if too many points
        npts = length(t)
        if npts > PLOT_MAX_POINTS
            step = max(1, Int(floor(npts / PLOT_MAX_POINTS)))
            idx = 1:step:npts
            t = t[idx]; S = S[idx]
        end

        Nval = getval(r, "N")
        lbl = Nval === nothing ? "entry_$i" : "N=$(Nval)"
        col = palette[mod1(i, length(palette))]
        plot!(p, t, S, label = lbl, color=col, linewidth=2.2)
    end

    xlims!(p, 0, maximum(map(r->maximum(Vector{Float64}(r.scaled_t)), results)))
    xlabel!(p, "Cycles (dimensionless)")

    out = replace(ARGS[1], ".jld2" => "_entropy_vs_N.pdf")
    savefig(p, out)
    println("Saved: $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
