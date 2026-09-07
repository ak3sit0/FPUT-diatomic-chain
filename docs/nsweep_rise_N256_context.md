# Config: nsweep_rise_N256

**Estado: ejecutado, resultado sin anotar.** Existen `configs/production/nsweep_rise_N256.toml`,
datos en `results/data/nsweep_rise_N256/` y figuras en
`results/figures/nsweep/entropy_Nsweep_N256.{png,pdf}`, pero este doc solo registra
la pregunta/diseño previo al run — falta anotar si Δκ=0.05/0.30 se mantuvieron planos
y si la subida en Δκ≈0.1–0.2 persiste a N=256 (criterios de éxito de la última sección).

## Pregunta física

¿Persiste para N>64 la subida parcial de entropía observada en N=64 a scaled_t ~1e3-1e5?

## Referencia (N=64)

Datos: `sweep_results_2026-07-02_periodic_case.jld2`, N=64, E_total=0.445, ε=0.00695

- Δκ=0.10: S/logN: 0.342 → 0.358 → 0.479 (scaled_t 1e3 → 1e4 → 6e4)
- Δκ=0.20: S/logN: 0.379 → 0.417 → 0.475
- Δκ=0.05, 0.30: planas en ~0.33 (controles, no suben)

**Interpretación:**
- Subida específica de Δκ ≈ 0.1–0.2
- E_opt/E ≈ 0 en todos los Δκ → la banda óptica no se llena
- Subida = redistribución DENTRO de la banda acústica, no transferencia interbanda
- Ningún Δκ alcanza equipartición (S/logN ≤ 0.48)
- Observable: altura e instante de la subida parcial, no un cruce de umbral alto

## Control del barrido en N

**ε FIJA** en 0.00695 (= 0.445/64, la del run de referencia)
- Implica E_total = ε·N en cada caso
- Alternativa (fijar E_total): ε caería como 1/N y los sistemas grandes serían progresivamente más lineales; la subida se debilitaría por eso y no por efecto de tamaño finito

**scaled_t_max FIJO**
- scaled_t = t·ω₂/2π con ω₂ ∝ 1/N ⇒ TMAX ∝ N
- Alternativa (TMAX fijo): N=256 apenas llegaría a scaled_t ≈ 3.9e4 y cortaría la subida por la mitad

**n_samples FIJO**
- nt constante ⇒ modal_E ∝ N (no ∝ N²)
- Evita explosión de memoria (modal_E sería ~100 MB por caso en N=256 si escalara como N²)

## Criterios de éxito

- Δκ=0.05 y 0.30 deben permanecer planos al crecer N
- Si aparece subida en ellos, el efecto no sería el de la referencia
