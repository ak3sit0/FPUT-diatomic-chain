# plot_entropy_comparison_N.jl
# Grafica entropía para distintos N en la misma figura para comparación

using JLD2, Plots, LaTeXStrings

include("../src/energy_analysis.jl")
using Main.EnergyAnalysis

# ============================================================================
# CONFIGURACIÓN 
# ============================================================================
DATA_DIR = "results/raw/alpha_sweep_N_periodic"
OUTPUT_DIR = "results/figures/entropy"
DELTA_SMOOTH = 0.6
NS = [32, 64, 128, 256]

mkpath(OUTPUT_DIR)

# ============================================================================
# CARGAR DATOS
# ============================================================================
println("📂 Buscando archivos .jld2 en $DATA_DIR...")

results_by_N = Dict{Int, Any}()

for N in NS
    # Buscar archivo con N específico
    files = filter(f -> contains(f, "N$N") && endswith(f, ".jld2"),
                   readdir(DATA_DIR))
    
    if isempty(files)
        println("  ⚠️  N=$N: No se encontró archivo")
        continue
    end
    
    filepath = joinpath(DATA_DIR, files[1])  # Tomar el primero
    println("  ✓ N=$N: Cargando $(basename(filepath))")
    
    try
        data = jldopen(filepath, "r")
        results = data["results"]
        
        # results es un vector de (param, delta, t, energies) tuples
        # Extraer las energías del primer resultado
        if !isempty(results)
            modal_E = results[1].modal_E
            scaled_t = results[1].scaled_t
            results_by_N[N] = (modal_E, scaled_t)
            println("    - Dimensiones: $(size(modal_E))")
        end
        
        close(data)
    catch e
        println("    ✗ Error al leer: $e")
    end
end

if isempty(results_by_N)
    println("\n❌ No se pudieron cargar datos. Verifica los archivos .jld2")
    exit(1)
end

# ============================================================================
# CALCULAR ENTROPÍA Y GRAFICAR
# ============================================================================
println("\n📊 Calculando entropías...")

p = plot(xlabel=L"t \, (\mathrm{ciclos})", ylabel=L"S \, (k_B)", 
         title=L"\alpha\text{-FPUT periodic BC, } \delta\kappa=0.1",
         legend=:bottomright, size=(900, 600), margin=5Plots.mm)

colors = [:red, :blue, :green, :orange]

for (idx, N) in enumerate(sort(collect(keys(results_by_N))))
    modal_E, scaled_t = results_by_N[N]
    entropy = compute_entropy(modal_E, DELTA_SMOOTH)
    
    # Plotear
    plot!(p, scaled_t, entropy; label=L"N = $N", linewidth=2.5, 
          color=colors[idx], alpha=0.8)
    
    println("  ✓ N=$N: S(0)=$(round(entropy[1]; digits=3)), " *
            "S(∞)≈$(round(entropy[end]; digits=3))")
end

# Exportar
output_path = joinpath(OUTPUT_DIR, "entropy_comparison_N.pdf")
savefig(p, output_path)
println("\n✅ Gráfico guardado en: $output_path")
