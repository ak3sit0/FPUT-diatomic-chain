# Diagnósticos Físicos — Notas para Futuro Trabajo

## Resonancia a Δκ = 0

**Estado:** Correcto, pero requiere vigilancia.

En `examples/resonance/plot_resonance_level_curves.jl`, cuando `Δκ = 0`, la cadena diatomic con κ₁ = κ₂ = 1 NO es una cadena monoatómica. Las relaciones de dispersión se acoplan, pero las curvas de nivel ω₋ + ω₋ = ω₊ quedan vacías (sin ceros). Esto está implementado correctamente en `compute_residual()` forzando `delta == 0.0 → fill(1.0, ...)` para evitar artefactos numéricos en el contour plotter.

**Atención:** Si se añade lógica de límite monoatómico explícito, verificar que no contradiga el tratamiento actual de Δκ → 0.

## Convención de Momentum y Wrapping

**Estado:** Verificado coherente en todo el código.

- `k3 = wrap_pi(-k1 - k2)` usado en `compute_gamma`, `resonance_matrix`, `plot_resonance_level_curves`
- `wrap_pi(k) = mod(k + π, 2π) - π` → mapea a [-π, π]
- Consistente con momentum conservation modular en lattice periódica

## Simetría de Rama en Coeficientes de Acoplamiento

**Estado:** Verificado tras Fix 2.

Bajo la parametrización de rama (s₁, s₂, s₃):
- Combinaciones simétricas (s₁ = s₂): aaa, aao, ooa, ooo → Γ(k₁,k₂) = Γ(k₂,k₁)
- Combinaciones asimétricas (s₁ ≠ s₂): aoa, oaa, aoo, oao → Γ(k₁,k₂) ≠ Γ(k₂,k₁)

El Fix 2 de orientación de matriz garantiza que esta simetría se refleja correctamente en los heatmaps de `plot_gamma_with_resonance.jl`.

## Energía e Entropia — Escala de Tiempo

**Referencia:** memory `entropy_rise_is_intra_acoustic` y `FPUT_cost_model`.

- Termalización observada es intraaústica, no transferencia interband
- Específica a Δκ ≈ 0.1–0.2 (ver `plot_entropy_Nsweep.jl`, `plot_thermalization_time.jl`)
- TMAX debe escalar ∝ N para resolver dinámicas resonantes a diferentes tamaños de sistema

## Doble formula de dispersión ω(k) — riesgo de divergencia silenciosa

**Estado:** Verificado (2026-08-16) — **no divergen hoy**; el riesgo es futuro.

Comparación numérica sobre k ∈ [-π,π] (1201 puntos) y Δκ ∈ {0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9}:

| Comparación | max \|Δω\| |
|---|---|
| `FPUTCoupling.omega_branch` vs `Dispersion.omega_ac`/`omega_op` | **0.0** (bit-idénticas) |
| forma compacta (`omega_compact_±`) vs forma absoluta | ~1e-15 |
| fórmula local de `plot_dispersion_relation.jl` vs ambas | ~1e-15 |

La tercera copia (la que vivía inline en `examples/dispersion/plot_dispersion_relation.jl` y en
`examples/resonance/plot_resonance_level_curves.jl`) **ya fue eliminada**: ambos scripts llaman
ahora a `src/dispersion.jl`. Quedan dos implementaciones (`Dispersion` y `FPUTCoupling.omega_branch`).

`src/dispersion.jl` (`omega_ac`, `omega_op`) y `src/fput_coupling.jl` (`omega_branch`) calculan la misma física — la dispersión ω(k) de la cadena diatómica — con fórmulas escritas independientemente en dos parametrizaciones distintas (κ₁,κ₂ absolutos vs Δκ compacto). No están unificadas: un cambio futuro en una convención (p.ej. signo de disc, definición de κ*) puede divergir silenciosamente de la otra sin que ningún test lo detecte, porque cada módulo se usa en scripts distintos y nada los compara entre sí.

**Sugerencia para más adelante:** un test de regresión que compare `omega_branch(k,kA,kB,·)` contra `omega_ac`/`omega_op` (o su forma compacta) en una grilla de k y Δκ, verificando equivalencia algebraica vía `sp.simplify` antes de confiar en que ambas parametrizaciones son intercambiables.

## Inconsistencia de esquema entre `compute_ensemble.jl` y `compute_ensemble_Nsweep.jl`

`compute_ensemble_Nsweep.jl` no calcula `T_therm_*` (tiempo de termalización), a diferencia de `compute_ensemble.jl`. Los scripts de plotting que consumen ambos formatos (`plot_ensemble_results.jl`, `plot_thermalization_time.jl`) asumen implícitamente que ese campo existe — si algún día se apunta `plot_thermalization_time.jl` a un resultado de `compute_ensemble_Nsweep.jl`, fallará o dará datos incompletos sin previo aviso. No es un bug de física en sí, pero afecta la interpretación de qué corridas son comparables entre sí.


## Dos definiciones distintas de "tiempo de termalización"

**Severidad: media** — afecta a la interpretación de las figuras, no a las simulaciones.

Conviven dos estimadores distintos de T_th bajo el mismo nombre:

| Dónde | Definición |
|---|---|
| `compute_ensemble.jl` → campos `T_therm_*` del JLD2 | umbral sobre la **entropía normalizada** por realización |
| `PlottingUtils.thermalization_time`, usado por `plot_ensemble_results.jl` y `plot_thermalization_time.jl` | primer t en que **E_optical** alcanza el 90% de su valor asintótico |

Consecuencia concreta: `plot_ensemble_results.jl` **ignora los campos `T_therm_*` que el JLD2 ya
trae** y recomputa T_th con la otra definición. Las dos cantidades no son comparables entre sí, y
la figura "Thermalization time vs Δκ" no muestra lo que `compute_ensemble.jl` registró.

Dado el hallazgo de memoria de que la subida de entropía es **intra-acústica** (E_opt/E ≈ 0 en
todos los Δκ para los casos de referencia), la definición basada en E_optical es sospechosa en ese
régimen: si la banda óptica nunca se llena, "el 90% del valor asintótico de E_opt" mide ruido.
Conviene decidir cuál es la definición canónica antes de usar estas figuras en el paper.
