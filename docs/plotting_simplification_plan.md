# Plan: adelgazar los scripts de plotting sin tocar la estética

Estado: **propuesto, no implementado**. Continuación de `plotting_backend_diagnostics.md`.

## Motivación

Los scripts de `scripts/plot/` son largos e incómodos de leer. Tras el refactor previo (paletas
centralizadas en `src/plotting_utils.jl`, dos scripts fusionados, duplicación de
`compute_entropy`/`thermalization_time` eliminada), lo que queda largo **no es lógica**: es
estética repetida línea a línea dentro de cada figura.

| Script | Líneas | Diagnóstico |
|---|---|---|
| `plot_ensemble_results.jl` | 230 | **una sola función de 190 líneas** — el problema real |
| `plot_heatmap_grid.jl` | 220 | ya son 9 funciones pequeñas; su estética vive en `PlotConfig`/TOML |
| `plot_ftmle.jl` | 140 | `fit_powerlaw` es lógica real, justificado |
| `plot_entropy_param_sweep.jl` | 135 | ya bien factorizado |
| `plot_entropy_paper.jl` | 98 | ya reducido |

Objetivo: máxima simplicidad **conservando el aspecto exacto** de las figuras. Todo el trabajo
salvo una excepción acordada debe salir byte-idéntico.

**Hallazgo que motiva la capa de tema:** los tres scripts de Plots.jl ya divergieron en silencio
en tipografía — `tickfont` 12 / 11 / 13 y `guidefont` 14 / 14 / 18 entre `plot_entropy_paper.jl`,
`plot_ftmle.jl` y `plot_entropy_param_sweep.jl`. Hoy nadie lo notaría. Centralizarlos no solo
acorta: pone las tres variantes juntas donde la divergencia se ve.

## Decisiones acordadas

1. **Tema primero, backend después.** Se mantienen Plots.jl y CairoMakie con dos capas de tema.
   La unificación (Fase 3) queda como paso siguiente sobre una base ya limpia, sin mezclar
   refactor estructural con cambio visual.
2. **Funciones nombradas**, no lista declarativa. Lectura lineal, sin indirección.
3. **Quitar el panel `E_opt/E_tot` repetido de la figura 4** de ensemble. Es copy-paste de la
   figura 3. La figura 4 pasa a un solo panel (`T_th`). *Es la única figura que cambia.*

## Cambios

### 1. Dos capas de tema (archivos nuevos)

Deliberadamente separadas para no forzar ambas dependencias en cada script: hoy
`plot_ensemble_results.jl` carga CairoMakie sin Plots y `plot_entropy_paper.jl` al revés.
Fundirlas en un solo archivo duplicaría el tiempo de arranque de todos.

**`src/theme_makie.jl`** (módulo `ThemeMakie`, `using CairoMakie`)

- `fput_theme()` → `Theme(...)` con `Axis` (titlesize, xlabelsize, ylabelsize), `Colorbar`
  (width, ticklabelsize) y `figure_padding`, para aplicar con `set_theme!`.
- `save_figure(fig, dir, name; px_per_unit=2)` — `save` + `println("Saved: …")`, hoy repetido
  4 veces en `plot_ensemble_results.jl`.
- `panel_row(fig, n; titles, ylabel)` — el esqueleto `Figure` + `GridLayout` + `Axis` por columna
  + `Label` inferior + `colgap!`, idéntico entre las figuras 1 y 2 de ensemble.
- **No** absorbe la estética de `plot_heatmap_grid.jl`: ese script ya la tiene externalizada en
  `PlotConfig` (`titlesize`, `ticksize`, `font`, `labelsize` vienen del TOML). Un tema global
  competiría con la config. Ese script solo adopta `save_figure`.

**`src/theme_plots.jl`** (módulo `ThemePlots`, `using Plots`)

- Tres presets nombrados que preservan los valores actuales *exactos* de cada script
  (`STYLE_PAPER`, `STYLE_SWEEP`, `STYLE_FTMLE`) y `apply_style!(preset)`. No se unifican los
  valores: eso cambiaría las figuras. Se ponen juntos para que la divergencia sea visible y se
  pueda decidir después.
- `log_time_plot(; ylabel, kwargs...)` — construye el `plot` base con eje x log y los ticks `10^k`
  en LaTeX (`xticks = (10.0 .^ (0:6), [L"10^{k}"])`), que hoy se reconstruyen a mano en
  `plot_entropy_paper.jl` y en `build_base_plot` de `plot_entropy_param_sweep.jl`. `size`, `grid`
  y márgenes quedan como argumentos explícitos porque **difieren** entre ambos (700×500 con grid
  punteado vs 800×600 sin grid).

### 2. `scripts/plot/plot_ensemble_results.jl` — romper el monolito

Es el grueso del trabajo. `plot_ensemble_diagnostics` (líneas 39–225) pasa a:

- `collect_stats(results)` → NamedTuple de vectores. Reemplaza cinco `push!` en paralelo seguidos
  de cinco reindexados con el mismo `sortperm` (líneas 60–87) por **un solo**
  `sort(results; by = r -> r.Delta)` y comprensiones. Reutiliza
  `PlottingUtils.thermalization_time`.
- `fig_modal_heatmaps(stats, cases)` — figura 1, sobre `panel_row`.
- `fig_entropy_bands(stats, cases)` — figura 2, sobre `panel_row`.
- `fig_global_summary(stats)` — figura 3, con los `delta_axis`/`trend!` que ya existen.
- `fig_thermalization(stats)` — figura 4, **ahora de un solo panel**.
- `print_summary_table(stats)` — la tabla final.
- `main` queda como cinco llamadas legibles.

Aplanar también el anidamiento `for → if haskey → if !isempty` (3 niveles) a guardas tempranas.

### 3. Los otros cuatro scripts

Cambios pequeños y mecánicos: sustituir su bloque de estilo por `apply_style!(STYLE_*)` o
`set_theme!(fput_theme())`, su `plot(...)` base por `log_time_plot(...)`, y sus `save`+`println`
por `save_figure`. No se reestructura nada más: ya están bien factorizados.

### 4. Documentación

- `plotting_backend_diagnostics.md`: registrar la divergencia tipográfica, y que la Fase 3 ahora
  arranca desde una base con la estética aislada en dos temas.
- `scripts/README.md`: mencionar los dos módulos de tema y cuándo usar cada uno.

## Verificación

Criterio: **byte-idéntico** (`cmp`), el mismo que ya funcionó en el refactor anterior para las
figuras de entropía.

1. Guardar copias de referencia de las figuras actuales antes de tocar nada.
2. Tras cada cambio, regenerar y `cmp` contra su referencia:
   - `plot_entropy_paper.jl` con los tres JLD2 de sweep (FBC, PBC, delta09) → 3 figuras.
   - `plot_entropy_param_sweep.jl` con `sweep_results_2026-07-02_periodic_case.jld2`.
   - `plot_ftmle.jl` con los tres JLD2 de `results/data/ftmle/`.
   - `plot_heatmap_grid.jl` con `sweep_results_2026-07-02_periodic_case.jld2` (el caso FBC de
     10⁶ columnas tarda mucho y consume ~5 GB; basta el periódico).
3. **`plot_ensemble_results.jl` no tiene datos en disco.** Sintetizar un JLD2 de humo en
   `results/data/_smoke/` con los campos que consume (`Delta`, `scaled_t`, `entropy_mean`,
   `entropy_std`, `modal_E_mean`, `E_acoustic_mean`, `E_optical_mean`), verificar las figuras 1–3
   byte-idénticas y **revisar la figura 4 a ojo** (cambia por diseño). Borrar el directorio de
   humo al terminar.
4. `julia --project=. test/runtests.jl` → 6/6 y 5/5.

Si una figura no sale byte-idéntica sin ser la figura 4 de ensemble, es un error del refactor:
corregirlo, no aceptarlo.

## Fuera de alcance

- Unificación de backend (Fase 3) — acordado como paso siguiente.
- `PlotConfig.outdir_xi` es un campo muerto (`plot_entropy_param_sweep.jl` documenta curvas ξ que
  no calcula).
- Las dos definiciones distintas de T_th (entropía normalizada en `compute_ensemble.jl` vs energía
  óptica en los scripts de figura) siguen pendientes de decisión física.
