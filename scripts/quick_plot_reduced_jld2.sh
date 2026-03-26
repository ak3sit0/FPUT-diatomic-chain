#!/bin/bash
set -euo pipefail

SRC="results/raw/alpha_optical_N64_periodic/results_springs_alpha_periodic_2026-03-06.jld2"
DST="results/raw/alpha_optical_N64_periodic/results_springs_alpha_periodic_2026-03-06_reduced.jld2"

julia --project=. - <<J
using JLD2
f = "${SRC}"
@info "loading" f
d = load(f)
results = d["results"]
keep_deltas = Set([0.01, 0.1, 0.2])
results_f = filter(r -> r.Delta in keep_deltas, results)
for i in 1:length(results_f)
    res = results_f[i]
    t = res.scaled_t
    e = res.modal_E
    idx = 1:10:length(t)
    results_f[i] = (param=res.param, Delta=res.Delta, scaled_t=t[idx], modal_E=e[:,idx])
end
@info "saving reduced file" "${DST}"
@save "${DST}" config=d["config"] results=results_f
J

cat > configs/plot_heatmap_periodic_optical_deltas_quick.toml <<EOF2
[data]
input_file = "${DST}"

[plot]
t_min          = 1e-3
t_max          = 1.0e18
time_subsample = 1
max_delta_cols = 3
color_scale    = "log"

[filter]
delta_values = [0.01, 0.1, 0.2]

[output]
heatmaps_dir = "results/figures/heatmaps"
EOF2

julia --project=. examples/plot_heatmap_grid.jl configs/plot_heatmap_periodic_optical_deltas_quick.toml
