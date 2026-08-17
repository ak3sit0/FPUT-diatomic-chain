# Examples — FPUT Diatomic Chain Simulations

Organización temática de scripts de análisis y visualización.

## Estructura

```
examples/
├── coupling/              # Coeficientes de acoplamiento Γ y scattering
│   ├── coupling_coefficients.jl     # Heatmap 3×3 de |Γ_σ1σ2σ3| (±curvas de resonancia)
│   ├── plot_aao_delta_sweep.jl      # Barrido en Δκ de |Γ_aao|
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
- **`src/dispersion.jl`** — Dispersión analítica (`omega_ac`, `omega_op`, formas compactas, derivadas)
- **`src/plotting_utils.jl`** — Paletas, estilos y helpers de entropía compartidos (`cyc`, `PALETTE_DELTA`, `PALETTE_COUPLING`, `entropy_series`, …)

## Cambios recientes

- Movidas funciones duplicadas (`compute_gamma`, `fetch_e`, etc.) a módulo `FPUTCoupling` en `src/fput_coupling.jl`
- Reorganización temática para claridad de propósito
- Ver `docs/plotting_backend_diagnostics.md` para problemas conocidos de backend y paleta
- **Fix orientación de matriz (verificado con `Contour.jl`):** `resonance_matrix` y el llenado
  de `Gamma` en `compute_gamma` usaban `D[j,i]`/`Gamma[idx][i2,i1]`, transpuesto respecto a la
  convención de Makie/Contour.jl (`Z[i,j] ↔ (x[i], y[j])`). El bug era invisible en las 4
  combinaciones de rama simétricas (`s1=s2`: aaa, aao, ooa, ooo) pero producía ejes k1↔k2
  intercambiados en las 4 asimétricas (aoa, oaa, aoo, oao).
  Corregido a `D[i,j]`/`Gamma[idx][i1,i2]`.
- **Optimización de `compute_gamma`:** `fetch_e` reemplazado por `fetch_e!` con buffer
  preasignado + `@view` en las columnas de eigenvectores del loop interno. Resultado idéntico
  numéricamente; ~125× menos allocations y ~1.8× más rápido (Nk=301, Ngrid=101).

## Notas de desarrollo

- Los scripts de `coupling/` muestran heatmaps 2D; no requieren visualización interactiva (CairoMakie, guardado a PDF/PNG)
- `plot_gamma_with_resonance.jl` se fusionó en `coupling_coefficients.jl`: era el mismo grid 3×3
  más los contornos de resonancia, ahora bajo la bandera `plot_gamma(...; resonance=true)`
- `plot_dispersion_relation.jl` y `plot_resonance_level_curves.jl` ya no reimplementan ω±(k):
  llaman a `src/dispersion.jl` (verificado idéntico a ~1e-15)
- Los scripts de `dispersion/` y `resonance/` pueden migrar a CairoMakie para uniformidad (actualmente Plots.jl)
- Todos los scripts incluyen guardias `if abspath(PROGRAM_FILE) == @__FILE__` para permitir `include()` sin efectos secundarios
