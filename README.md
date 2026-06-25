# Dynamics of Bose Systems: Resonances & Thermalization

Análisis numérico de la dinámica de sistemas de Bose mediante simulación de cadenas Fermi-Pasta-Ulam-Tsingou (FPUT) con énfasis en **resonancias de dispersión**, **relaciones de dispersión** y **dinámica de termalización**.

## Contenido del Proyecto

### Core Modules (`src/`)
- **`fput_core.jl`**: Definiciones físicas puras (masas, resortes, Hamiltonianos)
- **`fput_fast_runner.jl`**: Integración numérica con solvers de alta precisión (DifferentialEquations.jl)
- **`fput_analysis.jl`**: Análisis modal, entropía espectral, velocidades de dispersión
- **`config.jl`**: Gestión de configuración desde TOML

### Visualization Scripts (`examples/`)
- `plot_resonance_level_curves.jl` — Curvas de nivel de resonancia ω₋(k₁) + ω₋(k₂) = ω₊(k₃)
- `plot_dispersion_relation.jl` — Relaciones de dispersión ω±(k)
- `plot_gamma_with_resonance.jl` — Ancho de resonancia (gamma) en función de parámetros
- `compute_scattering_rate.jl` — Velocidades de dispersión (scattering rate)
- `plot_aao_delta_sweep.jl` — Parrón de barrido en parámetro delta
- `plot_entropy_N_timeseries.jl` — Evolución de entropía vs. tiempo
- `plot_thermalization_time.jl` — Análisis de tiempos de termalización
- `coupling_coefficients.jl` — Coeficientes de acoplamiento

### Simulation & Analysis (`scripts/`)
- `compute_trajectories.jl` — Ejecución de simulaciones FPUT
- `generate_hpc_jobs.jl` — Generación de scripts PBS para cluster
- `plot_entropy_param_sweep.jl` — Análisis paramétrico de entropía

## Quick Start

### Installation

```bash
git clone <repo-url>
cd Codigo
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

### Plotting Examples

```bash
# Curvas de resonancia
julia --project=. examples/plot_resonance_level_curves.jl

# Relación de dispersión
julia --project=. examples/plot_dispersion_relation.jl

# Análisis de velocidades de dispersión
julia --project=. examples/compute_scattering_rate.jl

# Termalización
julia --project=. examples/plot_thermalization_time.jl
```

## Directory Structure

```
├── configs/               # Archivos TOML para configuración
├── examples/              # Scripts de visualización y análisis
├── src/                   # Módulos core (simulación, análisis)
├── scripts/               # Scripts auxiliares y HPC
├── results/               # Salidas (figuras, datos)
│   ├── figures/           # Gráficas (PNG, PDF)
│   ├── raw/               # Datos simulados (.jld2)
│   └── logs/              # Logs de simulación
├── jobs/                  # PBS scripts generados
├── Project.toml           # Manifest de dependencias Julia
└── Manifest.toml
```

## Configuration

Los archivos TOML en `configs/` especifican parámetros físicos:

```toml
[physical]
N = 256                      # Número de partículas
DeltaK = 0.1                 # Parámetro de anisotropía
alpha = 0.5                  # Acoplamiento no-lineal

[output]
base_dir = "results/figures/resonance_level_curves"
```

## Dependencies

- `DifferentialEquations.jl` — Solvers ODE
- `Plots.jl` — Visualización
- `CairoMakie.jl` — Gráficos de alta calidad
- `TOML.jl` — Parseo de configuración
- `LaTeXStrings.jl` — Etiquetas LaTeX

## License

MIT
