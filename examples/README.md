# Examples — FPUT Diatomic Chain Simulations

Organización temática de scripts de análisis y visualización.

## Estructura

```
examples/
├── coupling/              # Coeficientes de acoplamiento Γ y scattering
│   ├── coupling_coefficients.jl     # Heatmap 3×3 de |Γ_σ1σ2σ3|
│   ├── plot_aao_delta_sweep.jl      # Barrido en Δκ de |Γ_aao|
│   ├── plot_gamma_with_resonance.jl # |Γ| con curvas de resonancia
│   └── compute_scattering_rate.jl   # R(Δκ) por co-area formula
│
├── dispersion/            # Relaciones de dispersión ω(k)
│   └── plot_dispersion_relation.jl  # Ramas acústica/óptica para varios Δκ
│
├── resonance/             # Condiciones de resonancia de 3 ondas
│   └── plot_resonance_level_curves.jl # Curvas de nivel ω₋+ω₋=ω₊
│
└── ensemble/              # Análisis de trayectorias y termalización
    ├── plot_entropy_Nsweep.jl       # S(t) para varios N
    └── plot_thermalization_time.jl  # T_th vs Δκ
```

## Uso

Cada script es autónomo y puede ejecutarse:

```bash
julia examples/coupling/coupling_coefficients.jl
julia examples/dispersion/plot_dispersion_relation.jl
julia examples/resonance/plot_resonance_level_curves.jl
julia examples/ensemble/plot_entropy_Nsweep.jl
```

Los scripts importan funciones centralizadas desde `src/`:
- **`src/fput_coupling.jl`** — `compute_gamma`, `compute_eigenvectors`, `resonance_matrix`, `omega_branch`
- **`src/dispersion.jl`** — Dispersión analítica (`omega_ac`, `omega_op`, derivadas)

## Cambios recientes

- Movidas funciones duplicadas (`compute_gamma`, `fetch_e`, etc.) a módulo `FPUTCoupling` en `src/fput_coupling.jl`
- Reorganización temática para claridad de propósito
- Ver `docs/plotting_backend_diagnostics.md` para problemas conocidos de backend y paleta

## Notas de desarrollo

- Los scripts de `coupling/` muestran heatmaps 2D; no requieren visualización interactiva (CairoMakie, guardado a PDF/PNG)
- Los scripts de `dispersion/` y `resonance/` pueden migrar a CairoMakie para uniformidad (actualmente Plots.jl)
- Todos los scripts incluyen guardias `if abspath(PROGRAM_FILE) == @__FILE__` para permitir `include()` sin efectos secundarios
