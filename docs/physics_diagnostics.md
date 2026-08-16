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
