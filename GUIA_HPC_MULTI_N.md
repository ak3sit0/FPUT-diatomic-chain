# Guía: Sistema HPC Multi-N para Ensambles FPUT

## 📋 Resumen Rápido

El sistema permite correr ensambles sobre **múltiples tamaños de sistema (N)** en paralelo en un clúster HPC con PBS/Torque. Cada valor de N se ejecuta como un job independiente con recursos escalados automáticamente.

```
Config maestra (N_values=[32,64,128,256])
           ↓
generate_hpc_jobs.jl
           ↓
jobs/ (4 TOML + 4 PBS)
           ↓
qsub jobs/ensemble_N*.pbs
           ↓
Resultados en results/data/ensemble_N_sweep_N*/ 
```

---

## 🚀 Flujo de Trabajo Paso a Paso

### **PASO 1: Revisar/Editar la Configuración Maestra**

Archivo: `configs/templates/ensemble_N_sweep.toml`

**Parámetros clave:**

```toml
[physics]
N_values = [32, 64, 128, 256]           # ← Qué tamaños de N correr
delta_values = [0.1, 0.15, ..., 0.5]   # ← Barrido de desorden
n_real = 50                             # ← Realizaciones por (N, delta)
param_values = [0.1]                    # ← Nonlineal (alpha=0.1)

[simulation]
TMAX = 1e6                              # ← Tiempo de simulación
T_block = 5e4                           # ← Tamaño de bloque temporal
DT = 0.05                               # ← Paso temporal
save_every = 100                        # ← Guardar cada N pasos
entropy_delta = 0.6                     # ← Parámetro entropía
```

**Cómo cambiar:**
- Añadir/quitar N: editar `N_values = [32, 64, 128]`
- Más realizaciones: cambiar `n_real = 100`
- Más tiempo: aumentar `TMAX = 1e7`

---

### **PASO 2: Generar los Job Files (Local, sin ejecutar)**

En tu máquina local (antes de ir al clúster):

```bash
cd /ruta/al/proyecto

# Generar todos los .toml y .pbs
julia --project=. scripts/generate_hpc_jobs.jl configs/templates/ensemble_N_sweep.toml
```

**Output esperado:**
```
=== Generador HPC: Ensamble Multi-N ===
Config: configs/templates/ensemble_N_sweep.toml
N_values: [32, 64, 128, 256]

✓ N=32   ppn=8, walltime=12:00:00
✓ N=64   ppn=8, walltime=12:00:00
✓ N=128  ppn=16, walltime=24:00:00
✓ N=256  ppn=16, walltime=48:00:00
```

**Archivos creados:**
```
jobs/
├── ensemble_N32.toml + ensemble_N32.pbs
├── ensemble_N64.toml + ensemble_N64.pbs
├── ensemble_N128.toml + ensemble_N128.pbs
└── ensemble_N256.toml + ensemble_N256.pbs
```

---

### **PASO 3: Revisar los Archivos Generados (Muy Importante)**

Antes de enviar al clúster, inspecciona:

```bash
# Ver estructura de un PBS
cat jobs/ensemble_N64.pbs

# Ver el TOML mutado (N=64 debe estar fijo)
cat jobs/ensemble_N64.toml | head -15
```

**Qué verificar en el PBS:**
- Nombre del job correcto: `#PBS -N ensemble_N64`
- Recursos: `ppn=8` o `ppn=16` según N
- Walltime: `12:00:00` o `24:00:00` o `48:00:00`
- Invocación: `julia --project=. -t <ppn> scripts/compute_ensemble.jl jobs/ensemble_N<N>.toml`

**Qué verificar en el TOML:**
- `N = 64` (debe ser número fijo, no array)
- NO debe tener `N_values`
- `base_dir = results/data/ensemble_N_sweep_N64` (con sufijo `_N<N>`)

---

### **PASO 4: Test Local (RECOMENDADO antes de HPC)**

Ejecutar un ensamble pequeño sin PBS (en tu máquina o login node):

```bash
# Test con N=32 (rápido, ~2-5 min depende del HW)
# Usa -t 4 si tienes 4 cores disponibles
julia --project=. -t 4 scripts/compute_ensemble.jl jobs/ensemble_N32.toml
```

**Qué significa que funcionó:**
```
[case 1] p=0.1  Δ=0.1  (init=band_ensemble, n_real=50)
  Banda 'acoustic': modos 2..3 (2 modos)
  Banda complementaria: modos 4..32 (29 modos)
  [r=1] S_final=...  E_ac=...  E_opt=...
  [r=2] S_final=...  E_ac=...  E_opt=...
  ...
[case 1] Finalizado. S̄=... ± ...
```

Y al final:
```
Guardado: results/data/ensemble_N_sweep_N32/ensemble_results_2026-06-24.jld2
```

---

### **PASO 5: Enviar a Clúster (PBS/Torque)**

**En el nodo de login del clúster:**

```bash
# Copiar los jobs/ a tu directorio de trabajo en el clúster
cd /tu/directorio/en/cluster
cp -r jobs/ .

# Enviar TODOS los jobs en paralelo
for job in jobs/ensemble_N*.pbs; do qsub $job; done
```

**Output esperado:**
```
12345.pbs-server
12346.pbs-server
12347.pbs-server
12348.pbs-server
```

Estos números son los **job IDs** — guárdalos si necesitas trackear.

---

### **PASO 6: Monitorear Progreso**

**En tiempo real:**
```bash
# Ver estado de tus jobs
qstat

# Ver logs (reemplaza 12345 con un job ID)
qstat -f 12345

# Seguir logs en directo
tail -f results/logs/ensemble_N*.log
```

**Output esperado:**
```
Job id                    Name             User            Time Use S Queue
------------------------- ---------------- --------------- -------- - -----
12345.pbs-server          ensemble_N32     joseangel       00:15:20 R default
12346.pbs-server          ensemble_N64     joseangel       00:08:45 R default
12347.pbs-server          ensemble_N128    joseangel       --        Q default
12348.pbs-server          ensemble_N256    joseangel       --        Q default
```

Estados:
- `Q` = En cola (Queued)
- `R` = Corriendo (Running)
- `E` = Error

---

### **PASO 7: Recolectar Resultados**

**Después de que todos terminen:**

```bash
# Ver outputs generados
ls -lh results/data/ensemble_N_sweep_N*/

# Esperado:
# results/data/ensemble_N_sweep_N32/ensemble_results_2026-06-24.jld2
# results/data/ensemble_N_sweep_N64/ensemble_results_2026-06-24.jld2
# results/data/ensemble_N_sweep_N128/ensemble_results_2026-06-24.jld2
# results/data/ensemble_N_sweep_N256/ensemble_results_2026-06-24.jld2

# Verificar tamaños (deben ser >100MB cada uno)
du -h results/data/ensemble_N_sweep_N*/*.jld2
```

---

### **PASO 8: Analizar Resultados**

**Plotear cada N individualmente:**

```bash
# Para N=32
julia --project=. examples/plot_thermalization_time.jl \
  results/data/ensemble_N_sweep_N32/ensemble_results_2026-06-24.jld2

# Para N=64
julia --project=. examples/plot_thermalization_time.jl \
  results/data/ensemble_N_sweep_N64/ensemble_results_2026-06-24.jld2
```

**Combinar múltiples N (requiere script adicional):**
- Cargar todos los JLD2 de diferentes N
- Plotear efectos de tamaño finito: cómo cambia T_th con N

---

---

## 📊 Sistema de Logging y Error Handling

### Cómo se capturan errores

Cada PBS script ahora:

1. **Genera un archivo log principal:**
   ```
   results/logs/ensemble_N64_<JOBID>.log
   ```
   Contiene TODO lo que Julia imprime (stdout + stderr mezclados)

2. **Genera marcadores de estado:**
   ```
   results/logs/ensemble_N64_<JOBID>.SUCCESS  ← si completó ok
   results/logs/ensemble_N64_<JOBID>.FAILED   ← si falló
   ```

3. **Captura el exit code:**
   El script falla si Julia falla (`set -e`), y reporta claramente

### Revisar logs después de ejecución

**Opción 1: Script automático (recomendado)**

```bash
julia --project=. scripts/check_hpc_status.jl
```

Output esperado:
```
================================================================================
HPC JOB STATUS REPORT
================================================================================

N       | Status               | Output     | Size
--------|----------------------|------------|------------------
32      | ✓ SUCCESS            | ✓ Exists   | 2.45 GB
64      | ✓ SUCCESS            | ✓ Exists   | 4.12 GB
128     | ✗ FAILED             | ✗ Missing  | —
256     | ⏳ RUNNING/QUEUED    | ✗ Missing  | —

================================================================================

[32] ✓ Completed successfully
     Log: ensemble_N32_12345.log
     Output: results/data/ensemble_N_sweep_N32/

[64] ✓ Completed successfully
     Log: ensemble_N64_12346.log
     Output: results/data/ensemble_N_sweep_N64/

[128] ✗ FAILED - Ver detalles:
     Failed marker: ensemble_N128_12347.FAILED
     Log: results/logs/ensemble_N128_12347.log

     --- LAST 30 LINES OF LOG ---
     Traceback error: segmentation fault in Julia thread...
     --- END OF LOG ---

[256] ⏳ Still running or in queue
     Log: results/logs/ensemble_N256_12348.log
     Ver progreso en tiempo real:
       tail -f results/logs/ensemble_N256_12348.log
```

**Opción 2: Ver logs manualmente**

```bash
# Ver último log de N=64
cat results/logs/ensemble_N64_*.log

# Ver solo últimas 50 líneas
tail -50 results/logs/ensemble_N64_*.log

# Seguir en tiempo real mientras corre
tail -f results/logs/ensemble_N64_*.log

# Buscar errores específicos
grep -i "error\|failed\|exception" results/logs/ensemble_N64_*.log
```

### Qué significa cada archivo marker

| Archivo | Significado | Acción |
|---------|------------|--------|
| `.SUCCESS` | Job completó sin errores | Analizar resultados |
| `.FAILED` | Job falló (ver log) | Revisar log, relanzar |
| Ninguno | Job aún no terminó o en cola | Esperar, ver con `qstat` |

### Si un job falla

**1. Ver el error en el log:**
```bash
tail -100 results/logs/ensemble_N128_12347.log
```

**2. Entender el error** (ejemplos comunes):

```
Error: OutOfMemoryError
→ Aumentar walltime o reducir n_real

StackOverflow / Segmentation fault
→ Bug en Julia o en el código FPUT

Killed by signal 9 (OOM killer)
→ Cluster mató el job por memoria, aumentar N_blocks

julia: not found
→ Julia no está en PATH, revisar ~/.bashrc en cluster
```

**3. Relanzar solo ese job:**
```bash
# Borrar archivos viejos (opcional)
rm results/logs/ensemble_N128_*.{SUCCESS,FAILED,log}
rm -rf results/data/ensemble_N_sweep_N128/

# Relanzar
qsub jobs/ensemble_N128.pbs
```

---

## 🛠️ Troubleshooting

### ❌ "No se generan los jobs"

```bash
# Verificar que la config maestra existe
ls -l configs/templates/ensemble_N_sweep.toml

# Verificar que tiene N_values
grep "N_values" configs/templates/ensemble_N_sweep.toml
```

### ❌ "Un job falla pero otros funcionan"

**Ventaja del diseño:** Los demás siguen.

Para relanzar solo N=128:
```bash
qsub jobs/ensemble_N128.pbs
```

### ❌ "¿Cuánto tiempo toma?"

Estimación (aproximado, depende del clúster):

| N | Cores | Walltime | Tiempo Real |
|---|-------|----------|-------------|
| 32 | 8 | 12h | 4-6h |
| 64 | 8 | 12h | 8-10h |
| 128 | 16 | 24h | 12-18h |
| 256 | 16 | 48h | 30-48h |

(Con n_real=50 realizaciones × 9 valores de delta)

### ❌ "Walltime muy poco/mucho"

Editar `configs/templates/ensemble_N_sweep.toml` y regenerar:
```bash
# Cambiar TMAX o n_real, luego:
julia --project=. scripts/generate_hpc_jobs.jl configs/templates/ensemble_N_sweep.toml
# Esto sobrescribe los jobs/ anteriores
```

---

## 📊 Estructura Final de Outputs

```
results/data/
└── ensemble_N_sweep_N32/
    ├── ensemble_results_2026-06-24.jld2  (N=32, 9 deltas, 50 real. c/uno)
    └── [contiene: Delta, scaled_t, E_optical_realizations, entropía, etc.]

└── ensemble_N_sweep_N64/
    ├── ensemble_results_2026-06-24.jld2  (N=64, misma estructura)
    ...

└── ensemble_N_sweep_N128/
└── ensemble_N_sweep_N256/
```

Cada JLD2 tiene la misma estructura; el N se distingue por el directorio.

---

## 💡 Tips Avanzados

### Cambiar solo algunos N

Editar manualmente `configs/templates/ensemble_N_sweep.toml`:
```toml
N_values = [64, 128]  # Solo correr estos dos
```

Luego regenerar.

### Cambiar recursos manualmente

Si los cores del clúster no son 8-16, editar en `scripts/generate_hpc_jobs.jl`:
```julia
function resources_for_N(N)
    if N <= 64
        return (ppn=12, walltime="08:00:00")  # ← cambiar aquí
    elseif N <= 128
        return (ppn=24, walltime="12:00:00")
    ...
end
```

Regenerar jobs.

### Dry-run (sin correr)

Ver qué pasaría sin ejecutar:
```bash
# Solo generar y ver estructura
julia --project=. scripts/generate_hpc_jobs.jl configs/templates/ensemble_N_sweep.toml

# Inspeccionar archivos
cat jobs/ensemble_N256.pbs
```

---

## ✅ Checklist Completo

- [ ] Revisar `configs/templates/ensemble_N_sweep.toml` (parámetros correctos)
- [ ] Generar: `julia scripts/generate_hpc_jobs.jl ...`
- [ ] Inspeccionar `jobs/*.pbs` y `jobs/*.toml`
- [ ] Test local: `julia -t 4 scripts/compute_ensemble.jl jobs/ensemble_N32.toml`
- [ ] Copiar `jobs/` al clúster
- [ ] Enviar: `for job in jobs/*.pbs; do qsub $job; done`
- [ ] Monitorear: `qstat`, `tail -f results/logs/`
- [ ] Analizar resultados cuando terminen

---

## 📞 Preguntas Frecuentes

**P: ¿Puedo cambiar walltime sin regenerar?**
No, está hardcodeado en `.pbs`. Editar config → regenerar.

**P: ¿Se pueden correr N en serie, no paralelo?**
Sí, mandar un job, esperar, mandar otro. Pero paralelo es lo normal.

**P: ¿Qué pasa si el clúster cae?**
Los jobs en cola se pierden. Los en ejecución pueden reencolarse.

**P: ¿Cómo combino resultados de múltiples N?**
Cargar cada JLD2, concatenar los DataFrames, plotear juntos. (Script a hacer si necesitas.)

---

¿Preguntas o algo que aclarar?
