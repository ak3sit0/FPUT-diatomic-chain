# Diagnóstico: Backend de Plotting y Paleta de Colores

## Problema 1: Coexistencia de backends sin justificación

**Estado actual:**
- **CairoMakie**: `coupling_coefficients.jl`, `plot_aao_delta_sweep.jl`, `plot_gamma_with_resonance.jl`, `compute_scattering_rate.jl`, `plot_entropy_Nsweep.jl`
- **Plots.jl**: `plot_dispersion_relation.jl`, `plot_resonance_level_curves.jl`, `plot_thermalization_time.jl`

**Raíz del problema:**
Históricamente, CairoMakie se eligió para gráficos 2D interactivos (heatmaps, contours de coeficientes). Plots.jl se usó para gráficos 1D (curvas de dispersión, líneas de resonancia). Sin embargo, ambos backends son capaces de producir ambos tipos de visualización, y la mezcla introduce inconsistencias de API y aumenta las dependencias del proyecto.

**Impacto:**
- Dos `savefig`/`save` con convenciones distintas
- Lógica de colormaps/paletas duplicada
- Los usuarios deben tener ambas librerías instaladas, aunque podrían usar una sola

## Problema 2: Paleta de colores inconsistente

**Síntoma:** Bug en `plot_dispersion_relation.jl:62`

```julia
col = colors[i]  # ← SIN mod1, causará BoundsError si i > 4
ls  = linestyles[mod1(i, length(linestyles))]
lw  = linewidths[mod1(i, length(linewidths))]
```

**Causa raíz:**
La mayoría de scripts usan `mod1()` para indexar arrays cíclicos de paleta, pero este script asume que `colors` tiene exactamente `length(Δκ_values_dispersion)=4` elementos. Cuando se añaden más valores de Δκ, el indexing falla silenciosamente hasta que i>4.

**Instancias de redefini ción de paleta:**
- `plot_gamma_with_resonance.jl`: hardcoded `:Blues`
- `plot_resonance_level_curves.jl`: paleta `[:darkblue, :steelblue, ...]` con `mod1`
- `plot_dispersion_relation.jl`: paleta `[:darkblue, :steelblue, ...]` sin `mod1`
- `plot_aao_delta_sweep.jl`: hardcoded `:Blues`
- `plot_entropy_Nsweep.jl`: paleta custom `[:royalblue, ...]`

No existe una única definición de "paleta estándar FPUT" reutilizable.

## Recomendaciones

### Corto plazo (semanas)
1. **Unificar a un solo backend**: CairoMakie es más moderno y permite tanto 2D como 1D sin compromiso de calidad. Mover `plot_resonance_level_curves.jl` a CairoMakie (requiere cambiar `contour()` → `contour!()` de Plots a `contours` de CairoMakie, que ya usa Contour.jl internamente).
2. **Centralizar paletas**: Crear `src/plotting_utils.jl` con:
   ```julia
   PALETTE_COUPLING = :Blues
   PALETTE_RESONANCE = [:darkblue, :steelblue, :cornflowerblue, :deepskyblue, :lightblue]
   PALETTE_DISPERSION = ... # mismo que resonance
   ```
3. **Fix inmediato**: Cambiar `colors[i]` → `colors[mod1(i, length(colors))]` en `plot_dispersion_relation.jl`.

### Mediano plazo (estándar futuro)
- Usar `Makie` (backend agnóstico) con un único `Makie.set_theme!()` global que unifica colormaps, tamaños de fuente, etc.
- Exportar figuras con `CairoMakie.activate!()` para PDF/PNG de alta calidad sin re-renderizar.

## Por qué no retraso esto más

La duplicación de paleta no es solo "estética". Con 8 combinaciones de rama en coeficientes acoplados + múltiples valores de Δκ en barridos, la indexación cíclica robusta (`mod1`) es **crítica** para evitar excepciones. El bug actual en `plot_dispersion_relation.jl` es el síntoma de que la lógica de indexing no está centralizada.
