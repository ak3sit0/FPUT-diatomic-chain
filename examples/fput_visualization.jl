"""
    fput_visualization.jl

Visualización literaria y funcional para las simulaciones FPUT.

Layers:
1. Data Loading (JLD2 -> Modal E)
2. Metric Computing (Entropy, Participation, xi)
3. Plotting Primitives (Heatmaps, Time-series)
4. Layout Orchestration (Master Plot)

Usage:
    julia --project=. examples/fput_visualization.jl results/raw/results_file.jld2
"""

using Plots, JLD2, Statistics, LaTeXStrings
include("../src/fput_core.jl");   using .FPUTCore
include("../src/fput_analysis.jl"); using .FPUTAnalysis

# ── Helpers ──
apply_style() = default(titlefont=font(14), guidefont=font(12), tickfont=font(10), legendfont=font(10), framestyle=:box)

# ── Layer 3: Plotting Primitives ──

function plot_energy_heatmap(modal_E, time; scale=:log10)
    data = (scale == :log10) ? log10.(modal_E .+ 1e-15) : modal_E
    h = heatmap(time, 1:size(modal_E, 1), data,
                xlabel=L"t", ylabel="Mode index",
                title="Modal Energy Distribution",
                color=:viridis, clims=(-12, 0))
    return h
end

function plot_entropy_series(time, entropy)
    p = plot(time, entropy, xscale=:log10,
             xlabel=L"t", ylabel=L"S(t)",
             title="Spectral Entropy Decay",
             lw=2, color=:crimson, label="S(t)")
    return p
end

function plot_resonance_level_curves(sys_params)
    # Accept either a saved Dict or a SystemParams struct
    if isa(sys_params, AbstractDict)
        sp = FPUTCore.SystemParams(Int(sys_params["N"]), Float64(get(sys_params, "DeltaK", 0.0)), Float64(get(sys_params, "DeltaM", 0.0)), Float64(get(sys_params, "alpha", 0.0)), Float64(get(sys_params, "beta", 0.0)), Symbol(sys_params["boundary"]))
    else
        sp = sys_params
    end

    # Recreate the dispersion relation omega_k vs mode k
    # Plotting resonance condition roughly Omega_i + Omega_j = Omega_k
    kvec, mvec = FPUTCore.make_system(sp)
    freq, _ = FPUTCore.find_normal_modes(kvec, mvec, sp.boundary)
    N = sp.N
    
    # Delta omega = |w_k - w_i - w_j|
    Z = zeros(N, N)
    target_mode = sp.boundary == :fixed ? 1 : 2 # Basic acoustic interaction assumption
    for i in 1:N
        for j in 1:N
            Z[i, j] = abs(freq[target_mode] - freq[i] - freq[j])
        end
    end
    
    p = contourf(1:N, 1:N, 1 ./ (Z .+ 1e-3),
                 title="Resonance Interaction (w_k - w_i - w_j)",
                 xlabel="Mode i", ylabel="Mode j",
                 colormap=:magma, levels=20)
    return p
end

# ── Layer 4: Visualization Narrative ──

function visualize_fput_results(input_path::String)
    println("--- Generating Literate Visualizations ---")
    data = load(input_path)
    modal_E = data["modal_E"]
    entropy = data["entropy"]
    time = data["time"]
    sys_params = data["config"] # Extract saved physical parameters
    
    apply_style()
    
    # Prepare individual paths
    base_out = replace(input_path, ".jld2" => "")
    
    # Generate Heatmap
    p1 = plot_energy_heatmap(modal_E, time)
    savefig(p1, "$(base_out)_heatmap.pdf")
    println("Saved: $(base_out)_heatmap.pdf")
    
    # Generate Entropy Plot
    p2 = plot_entropy_series(time, entropy)
    savefig(p2, "$(base_out)_entropy.pdf")
    println("Saved: $(base_out)_entropy.pdf")
    
    # Generate Resonance Level Curves
    p3 = plot_resonance_level_curves(sys_params)
    savefig(p3, "$(base_out)_resonance.pdf")
    println("Saved: $(base_out)_resonance.pdf")
    
    return p1, p2, p3
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        println("Uso: julia examples/fput_visualization.jl <results.jld2>")
    else
        visualize_fput_results(ARGS[1])
    end
end
