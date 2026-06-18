# Plan de cambios: ensamble de fases y barrido en Δκ

## Objetivo
Pasar de una condición inicial determinista (modo puro, `v=0`) a un **ensamble pequeño de fases aleatorias por banda**, para (1) probar que el crossover en Δκ es robusto a la condición inicial y no un artefacto de la fase nula, y (2) medir la transferencia interbanda en su forma estándar. Alinea el método con Pezzi/Onorato (diatómico) sin reescribir la física.

## Decisiones de diseño (cerradas en la discusión)
- **Promedio**: ensamble sobre fases (~8 realizaciones), **no** temporal. La ventana deslizante `μ=2/3` actual se conserva solo como suavizado dentro del plateau.
- **Rango físico**: Δκ ∈ {0.1, 0.15, 0.2, 0.3, 0.4, 0.5}. Excluir Δκ < 0.1 (band folding) y notar que el canal de 3 ondas muere en Δκ ≈ 0.5 (no 0.6: la condición de borde de banda da 3/5 pero la conservación de momento lo cierra antes).
- **Inicialización**: por banda selectiva (acústica u óptica), en **modos normales reales** (sin ángulo-acción).

---

## Cambio 1 — Inicialización de banda selectiva con fases (núcleo)
Reemplazar `q_cur .= amplitude .* U[:, target_mode]` por una función que:

1. Recibe `branch ∈ {:acoustic, :optical}`, una banda de índices de modo `k_band`, energía `ε`, y una semilla.
2. Fija el espectro inicial `E_k` (plano sobre `k_band` de la rama elegida, cero en la otra), normalizado a `Σ E_k = N·ε`.
3. Para cada modo real `j` de la banda, sortea `φ_j ~ U[0,2π)` y reparte:
   - `Q_j = sqrt(2 E_j)/ω_j · cos(φ_j)`
   - `P_j = -sqrt(2 E_j) · sin(φ_j)`
   (conserva `E_j` exactamente, independiente de φ).
4. Transforma modal→física con la matriz de autovectores ya disponible (`U` de `eigen`): `q = U·Q`, `p = U·P` (con la normalización de masa que ya usa el código).

**Por qué así**: la fase aleatoria es solo la rotación cos/sin del reparto Q/P dentro de cada oscilador modal; no requiere reescribir el Hamiltoniano en I-φ.

## Cambio 2 — Test de validación (barato, obligatorio)
Tras construir el estado inicial, proyectarlo de vuelta a la base modal y verificar:
- `E_{otra banda}(0) ≈ 0` (precisión de máquina) → confirma que la transformación inversa no tiene fuga.
- `Σ E_k(0) ≈ N·ε` → confirma normalización.
Dejar como `@assert` en la inicialización.

## Cambio 3 — Bucle de ensamble
Envolver la corrida en un bucle de `n_real = 8` semillas:
- Cada realización: misma `(Δκ, branch, ε, k_band)`, distinta semilla de fases.
- Acumular `S̄(t)` y energías modales por realización.
- Reportar **media y dispersión** entre realizaciones (la dispersión es el diagnóstico sticky vs resonante: estrecha → resonancia; cola larga → stickiness).

## Cambio 4 — Barrido en Δκ y heatmaps promediados
- Bucle externo sobre Δκ del rango físico.
- Para cada Δκ: ensamble (Cambio 3), guardar heatmap de energía modal vs tiempo promediado sobre el ensamble.
- Régimen atrapado (Δκ ≳ 0.5): **no** esperar termalización; medir saturación del plateau sobre la cola temporal.

## Cambio 5 — Lyapunov de tiempo finito (respaldo, distingue sticky vs resonante)
Sobre 2–3 realizaciones representativas por Δκ:
- Integrar trayectoria de referencia + copia desplazada `d0 ~ 1e-8`, renormalizar periódicamente (Benettin).
- `λ_ft(t) = (1/t) Σ log(d_i/d0)`.
- Firma: `λ_ft` decae ~1/t durante el plateau (regular) y salta al escapar (sticky). Correlacionar el salto con la subida de `S̄`.
- Hooks ahora; implementación completa solo si se confirma cola larga en Cambio 3.

---

## Orden de ejecución
1. Cambio 1 + 2 (inicialización + validación) — base de todo.
2. Cambio 3 + 4 (ensamble + barrido) — produce el resultado principal.
3. Cambio 5 (Lyapunov) — solo si la dispersión del ensamble sugiere stickiness.

## Qué NO cambia
- Integrador KahanLi8, `solve_fput`, cálculo de energías modales y entropía espectral: intactos.
- La inicialización determinista actual se conserva como caso `n_real=1, φ=0` para reproducibilidad.

## Predicciones contra las que contrastar
- Δκ > 0.5: atrapamiento (sin canal de 3 ondas).
- Δκ ∈ (0.1, 0.5): transferencia más rápida conforme Δκ → 0.1 (dominado por tamaño del manifold).
- Tiempo de transferencia ∝ 1/R_full(Δκ) calculado previamente.
