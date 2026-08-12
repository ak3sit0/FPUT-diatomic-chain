# nsweep_rise: ¿la subida de entropía persiste para N > 64?

## Lanzar

```bash
qsub jobs/nsweep_N64.pbs
qsub jobs/nsweep_N128.pbs
qsub jobs/nsweep_N256.pbs
```

Cada job corre los 4 Δκ (0.05, 0.1, 0.2, 0.3) en 4 hilos.

| job | walltime | salida |
|---|---|---|
| nsweep_N64 | 01:41 | `results/data/nsweep_rise_N64/` |
| nsweep_N128 | 06:47 | `results/data/nsweep_rise_N128/` |
| nsweep_N256 | 27:10 | `results/data/nsweep_rise_N256/` |

En local, los tres a la vez:

```bash
julia --project=. -t 12 scripts/compute_ensemble_Nsweep.jl configs/cases/nsweep_rise_N256.toml
```

## Regenerar los jobs

Si cambias el config (`configs/cases/nsweep_rise_N256.toml`):

```bash
julia --project=. scripts/generate_pbs_nsweep.jl configs/cases/nsweep_rise_N256.toml
```

El walltime sale del modelo de coste, no a ojo. La constante por defecto (67 ns/paso/sitio)
es la del portátil; para ajustarla al clúster:

```bash
julia --project=. scripts/calibrate_cost.jl        # → C = XX
julia --project=. scripts/generate_pbs_nsweep.jl <config> --ns-per-step-site=XX
```

## Qué mirar en los resultados

**Usar `exp(S)`, no `S/log N`.** `exp(S)` = número efectivo de modos excitados. Si el
paquete se ensancha hasta ~7 modos en toda N, `S/log N` baja de 0.48 (N=64) a 0.36 (N=256)
sólo por el denominador, y parecería que la subida no persiste siendo la física la misma.

Referencia N=64: `exp(S)` va de ~4 a ~7 en Δκ=0.1 y 0.2; se queda en ~4 en Δκ=0.05 y 0.3.
Toda la energía vive en los modos 2–11 y `E_opt/E ≈ 0` (la banda óptica no se llena).

| resultado | interpretación |
|---|---|
| `exp(S)` → ~7 en toda N | lo fija el índice de modo: la subida persiste igual |
| `exp(S)` ∝ N (7 → 14 → 28) | lo fija una fracción de la zona de Brillouin |
| `exp(S)` → ~4 | era un efecto de tamaño finito de N=64 |

Δκ=0.05 y 0.3 son controles: deben seguir planos.

## Notas

- ε = 0.00695 fija para toda N (= 0.445/64, la del run de referencia).
- `scaled_t_max = 1e6` ⇒ TMAX ∝ N. Con TMAX fijo, N=256 se quedaría en scaled_t≈3.9e4 y
  cortaría la subida por la mitad.
- Una trayectoria por caso, sin promedio de ensamble: `exp(S)` fluctúa. Si la comparación
  entre N sale marginal, hay que promediar sobre varias condiciones iniciales perturbadas.
- Coste ∝ N² a scaled_t fijo. N=512 serían ~25 h por Δκ.
