# REFACTORED: Using SimulationRunner for energy-based initialization test

include("../src/parameters.jl");       using Main.Parameters
include("../src/dynamical_matrix.jl"); using Main.DynamicalMatrix
include("../src/simulation_runner.jl"); using Main.SimulationRunner
include("../src/energy_analysis.jl");  using Main.EnergyAnalysis
include("../src/plotting.jl");         using Main.Plotting 

using LinearAlgebra, Plots 

println("Initializing simulation via SimulationRunner...")

# 1. Parámetros de la prueba
N = Parameters.N
target_energy = Parameters.initial_energy # Asegúrate de que existe en src/parameters.jl
init_mode = 1 # Excitar el primer modo

# 2. Ejecutar simulación y análisis
# SimulationRunner ahora calcula internamente: amplitude = sqrt(2*E) / omega
modal_E, scaled_t, freqs, v_final = SimulationRunner.run_and_analyze(
    N = N,
    TMAX = Parameters.TMAX,
    DT = Parameters.DT,
    alpha = Parameters.alpha, 
    beta = Parameters.beta,
    boundary = :fixed,
    init_mode = init_mode,
    energy = target_energy
)

println("Simulation complete. Target Energy: ", target_energy)
println("Energy at t=0 (Mode $init_mode): ", modal_E[init_mode, 1])

# 3. Verificación de energía total
# Sumamos las energías de todos los modos en cada instante de tiempo
total_E = sum(modal_E, dims=1)[:]

# 4. Generar Plots usando el módulo Plotting
println("Plotting results...")

# Plot de energía total (debe ser constante si alpha/beta son pequeños)
p1 = Plotting.plot_total_energy(scaled_t, total_E)
title!(p1, "Total Energy (Target: $target_energy)")

# Plot de evolución de energías modales
p2 = Plotting.plot_modal_energies(scaled_t, modal_E)
title!(p2, "Modal Energy Evolution (Init Mode: $init_mode)")

# Mostrar o guardar
#display(p1)
#display(p2)

savefig(p1, "results/figures/test_total_energy_initialization.png")
savefig(p2, "results/figures/test_modal_energy_initialization.png")