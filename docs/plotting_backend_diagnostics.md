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

**Instancias de redefinición de paleta:**
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

## Extensión: scripts/ (plots del paper y ensembles)

Mismo problema de fondo, con más superficie:

- **Split de backend**: `plot_ensemble_results.jl` y `plot_heatmap_grid.jl` usan CairoMakie; `plot_entropy_paper.jl`, `plot_entropy_param_sweep.jl`, `plot_entropy_pbc_complete.jl`, `plot_ftmle.jl` usan Plots.jl. Son figuras del **mismo paper** producidas por dos stacks distintos con manejo de fuente/DPI/márgenes diferente — riesgo real de inconsistencia visual entre paneles de una misma publicación.
- **Tres sistemas de color independientes** sin fuente única: los hex arrays `SPEC_FBC`/`SPEC_PBC`/`SPEC_PBC_COMPLETE` (`plot_entropy_paper.jl`, `plot_entropy_pbc_complete.jl`), el gradient `BLUES_GRADIENT` duplicado literalmente en `plot_ensemble_results.jl` y `plot_heatmap_grid.jl`, y el `PALETTE`/`LINESTYLES` de `plot_ftmle.jl`.
- Colores hex como `"#0D4CB3"` repetidos ~8 veces como literales en `plot_ensemble_results.jl` en vez de constantes nombradas — riesgo de que diverjan silenciosamente de `BLUES_GRADIENT`.
- Punto positivo: el patrón `colors[i]` sin `mod1` (el bug de `plot_dispersion_relation.jl`) **no** aparece en `scripts/` — `plot_entropy_param_sweep.jl` y `plot_ftmle.jl` ya usan `mod1` correctamente.

La recomendación de "unificar a un solo backend + centralizar paletas en `src/plotting_utils.jl`" aplica igual aquí, y cubriría examples/ y scripts/ a la vez.

## Estado: resuelto (paletas y duplicación)

Implementado en `src/plotting_utils.jl` — fuente única de paletas, estilos y helpers de figura:

| Antes | Ahora |
|---|---|
| `BLUES_GRADIENT` duplicado en `plot_ensemble_results.jl` / `plot_heatmap_grid.jl` | `BLUES_STOPS` |
| `SPEC_FBC`/`SPEC_PBC`/`SPEC_PBC_COMPLETE` en dos scripts | `SPEC_*` en `PlottingUtils` |
| `PALETTE`/`LINESTYLES` locales de `plot_ftmle.jl` | `FTMLE_PALETTE`/`FTMLE_LINESTYLES` |
| `:Blues`/`:orangered` hardcoded en 3 scripts de coupling | `PALETTE_COUPLING`/`RESONANCE_LINE` |
| `~8 literales "#0D4CB3"` en `plot_ensemble_results.jl` | `BLUE_DARK`/`BLUE_MID`/`GRAY_GUIDE` |
| paleta rainbow local en `plot_entropy_Nsweep.jl` | `PALETTE_CATEGORICAL` (mismos colores) |
| `palette[i]` / `colors[mod1(...)]` mezclados | `cyc(palette, i)` en todos los scripts |

**Deduplicación de lógica:**

- `compute_entropy` (duplicado byte-a-byte entre `plot_entropy_paper.jl` y
  `plot_entropy_pbc_complete.jl`) → `PlottingUtils.entropy_series`, que delega en
  `FPUTAnalysis.spectral_entropy`. Se verificó equivalencia numérica antes de sustituir:
  `max |Δ| = 6.6e-14`, `isapprox(rtol=1e-10) == true`. Las figuras regeneradas salen
  **idénticas byte a byte** a las anteriores.
- `prepare_ts` / `logdownsample` / `load_results_by_delta` → `PlottingUtils`.
- `plot_entropy_pbc_complete.jl` **eliminado**: era `plot_entropy_paper.jl` con otra spec y un
  loader multi-archivo. Ahora `plot_entropy_paper.jl <fbc> <pbc> [<pbc_extra>…]` produce las
  tres figuras.
- `plot_gamma_with_resonance.jl` **eliminado**: era `coupling_coefficients.jl` más contornos.
  Ahora `plot_gamma(kA, kB, α; resonance=true)`.
- `plot_dispersion_relation.jl` y `plot_resonance_level_curves.jl` reimplementaban ω±(k);
  ahora llaman a `src/dispersion.jl` (verificado idéntico a ~1e-15 sobre Δκ ∈ [0, 0.9]).

**Bugs reales encontrados al regenerar figuras** (no eran sólo estética):

1. `plot_heatmap_grid.jl` y `plot_entropy_param_sweep.jl` **crasheaban** con todos los JLD2
   existentes. Los resultados guardan la ruta del config (`configs/cases/…_production.toml`),
   que dejó de existir tras reorganizar `configs/`; el código sólo convertía la cadena a
   `ExperimentConfig` `if isfile(config)`, así que se quedaba `String` y reventaba en
   `config.nonlinear`. Resuelto con `Config.resolve_config_path` / `Config.as_experiment_config`
   (busca por basename en `configs/`, tolera el sufijo `_production` retirado, y degrada a un
   config por defecto con `@warn` en vez de lanzar).
2. `decade_ticks` (antes inline en `make_panel_axis`): la rama de respaldo construía etiquetas
   con `string(Int(x))` sobre potencias de diez que pueden ser < 1 → `InexactError`. Corregido.

## Rendimiento de `plot_heatmap_grid.jl`

`prepare_panel_data` troceaba la matriz modal en dos pasos
(`z[:, ventana][:, ::subsample]`), copiando la matriz N×nt completa en cada uno, más una
conversión `Float64.(modal_E)` previa. En una corrida de producción (64 × 10⁶ por panel, 6 paneles)
eso son ~0.5 GB por copia y ~2 GB vivos de forma transitoria por panel. Ahora se resuelven los
índices primero y se corta **una sola vez**; salida **bit-idéntica** (verificado con `cmp`).

**Intento descartado — dejar constancia para no repetirlo:** submuestrear el eje temporal a ~una
columna por píxel (2·cell_w·px_per_unit) parecía gratis, con el argumento de que el rasterizador
colapsa igual ~1700 columnas por píxel. **Es falso.** Medido contra la figura a resolución completa:

| Reducción | píxeles que difieren | >8/255 | RMS |
|---|---|---|---|
| escoger 1 por bin (log) | 15.8% | — | 4.61 |
| promediar por bin (log) | 18.5% | 6.5% | 4.78 |

Makie antialiasa las columnas sobrantes, así que ninguna reducción a una columna por píxel
reproduce la figura. Además el tiempo apenas bajó (46s → 42s en el caso de 2·10⁵ puntos): el coste
está en el arranque de Julia, la carga del JLD2 y el rasterizado, no en el número de columnas.
Se revirtió; queda sólo un `@warn` cuando un panel supera 2·10⁵ columnas, señalando `time_subsample`
como la vía consciente de bajar fidelidad a cambio de velocidad.

## Pendiente

- **Fase 3 — Unificación de backend**: siguen conviviendo Plots.jl (`plot_entropy_paper.jl`,
  `plot_entropy_param_sweep.jl`, `plot_ftmle.jl`, `plot_dispersion_relation.jl`,
  `plot_resonance_level_curves.jl`, `plot_thermalization_time.jl`) y CairoMakie
  (`plot_ensemble_results.jl`, `plot_heatmap_grid.jl`, `examples/coupling/*`). Ahora que las
  paletas están centralizadas y son agnósticas de backend, migrar es sustancialmente más barato,
  pero sigue siendo reescritura real de recipes (`contour`, `plot!`) y exige revisar visualmente
  cada figura regenerada.
- **Cambio de tono deliberado**: `plot_dispersion_relation.jl` usaba
  `[:darkblue, :steelblue, :royalblue, :cornflowerblue]`; ahora usa `PALETTE_DELTA`
  (`…, :cornflowerblue, :deepskyblue, …`), que es la que ya usaba `plot_resonance_level_curves.jl`.
  Las curvas 3 y 4 cambian de tono a cambio de que las dos figuras del paper sean consistentes.
  Revertir es cambiar una constante si se prefiere la paleta antigua.
- **`PlotConfig.outdir_xi`** es un campo muerto: `plot_entropy_param_sweep.jl` documenta curvas de
  localización ξ pero no las calcula.
