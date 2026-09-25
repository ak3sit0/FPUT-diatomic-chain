"""
    dimer_approach.jl

Approach to Δκ = 1 with the same initial condition as `single_dimer.jl`: only the
relative coordinate of dimer 1 excited (E = 0.445, zero centre of mass, every other
mass at rest). Same physical time for every Δκ, 10⁵ periods of ω = 2 (the Δκ = 1
dimer frequency), so the curves share a clock.

S̄(t) is computed in two bases, because they stop agreeing as soon as Δκ < 1:

- normal modes (the pipeline's S̄): for Δκ < 1 they are delocalized Bloch waves,
  so a single excited dimer already projects onto the whole optical band at t = 0;
- dimers (relative + centre-of-mass energy of each dimer, local in space): measures
  how far the energy has spread along the chain. The κ₂ bonds between dimers are
  not assigned to either, a (1−Δκ)-small share of the energy.

Initial condition (5th argument):
- `dimer` (default): the single-dimer excitation above.
- `mode2`: the production initial condition — all the energy in mode 2, the lowest
  acoustic mode. Undefined at Δκ = 1, where ω₂ = 0 (build_case refuses it too).
- `mode1`: the same with mode 1, the production choice for FBC.

Boundary (6th argument): `periodic` (default) or `fixed`. FBC outputs get `_fbc`.

Usage:
  julia --project=. dimer_limit/dimer_approach.jl [N=64] [E=0.445] [cycles=1e5] [deltas=0.9,0.95,0.99,1] [ic=dimer|mode2|mode1] [bc=periodic|fixed]

  julia --project=. dimer_limit/dimer_approach.jl --replot [--norm] dimer_limit/data/dimer_approach_*.jld2
  julia --project=. dimer_limit/dimer_approach.jl --merge a.jld2 b.jld2 …    # combine sweeps, plot S̄/ln N

Outputs are tagged with the Δκ list, so different sets do not overwrite each other.
"""

using LinearAlgebra, JLD2, Printf, LaTeXStrings, Plots
const ROOT = joinpath(@__DIR__, "..")
include(joinpath(ROOT, "src/fput_core.jl"));        using .FPUTCore
include(joinpath(ROOT, "src/fput_fast_runner.jl")); using .FPUTFastRunner
include(joinpath(ROOT, "src/fput_analysis.jl"));    using .FPUTAnalysis
include(joinpath(ROOT, "src/plotting_utils.jl"));   using .PlottingUtils
include(joinpath(ROOT, "src/plot_style.jl"));       using .PlotStyle

arg(i, d) = length(ARGS) >= i ? parse(Float64, ARGS[i]) : d

const DEFAULT_DELTAS = [0.9, 0.95, 0.99, 1.0]
const α, β, DT = 0.1, 0.0, 0.05
const ω_CLOCK = 2.0                   # Δκ = 1 dimer frequency, common clock for all runs
const NSAVE = 50_000                  # uniform: spectral_entropy's window is index-based

"""Total Hamiltonian (unit masses). FBC: N+1 bonds, the walls are fixed at q = 0."""
function hamiltonian(q, v, k, bc)
    d = bc == :fixed ? diff([0.0; q; 0.0]) : [q[mod1(i+1, length(q))] - q[i] for i in eachindex(q)]
    0.5 * sum(abs2, v) + sum(bond_potential.(k, d, α, β))
end

"""
Dimer-coordinate energies, N × nt (κ₁ = strong bond).

PBC: dimers are masses (2j−1, 2j); odd rows relative, even rows centre of mass.
FBC: bonds are wall–1 (κ₁), 1–2 (κ₂), 2–3 (κ₁), …, N–wall (κ₁), so the dimers are
(2,3), …, (N−2, N−1) and masses 1 and N are each tied to a wall by κ₁. Row 1 and
row N hold those two wall oscillators; rows 2…N−1 alternate relative / COM.
"""
function dimer_energies(Q, V, κ₁, bc)
    I = bc == :fixed ? (2:2:size(Q, 1)-2) : (1:2:size(Q, 1))     # first mass of each dimer
    r, ṙ = Q[I .+ 1, :] .- Q[I, :], V[I .+ 1, :] .- V[I, :]
    Q̇    = (V[I .+ 1, :] .+ V[I, :]) ./ 2
    E = similar(Q)
    E[I, :]      .= ṙ .^ 2 ./ 4 .+ κ₁ .* r .^ 2 ./ 2
    E[I .+ 1, :] .= Q̇ .^ 2
    if bc == :fixed
        for w in (1, size(Q, 1))
            E[w, :] .= V[w, :] .^ 2 ./ 2 .+ κ₁ .* Q[w, :] .^ 2 ./ 2
        end
    end
    E
end

"""Projected mode for the modal initial condition: PBC mode 2, FBC mode 1 (project convention)."""
const IC_MODE = Dict("mode2" => 2, "mode1" => 1)

function run(dk, N, E_in, T, ic, bc)
    sp   = SystemParams(N, dk, 0.0, α, β, bc)
    k, m = make_system(sp)
    freq, Vm = find_normal_modes(k, m, bc)

    q0, v0 = zeros(N), zeros(N)
    if ic == "dimer"
        r0 = sqrt(2 * E_in / (1 + dk))          # harmonic amplitude on the κ₁ bond
        i  = bc == :fixed ? 2 : 1               # FBC: dimer 1 = masses (2, 3)
        q0[i], q0[i+1] = -r0 / 2, r0 / 2
    elseif haskey(IC_MODE, ic)
        n = IC_MODE[ic]
        ω = ref_frequency(freq, n)              # errors at Δκ = 1 (acoustic band at 0)
        q0 .= sqrt(2 * E_in) / ω .* (Diagonal(1 ./ sqrt.(m)) * Vm)[:, n]
    else
        error("ic must be \"dimer\", \"mode1\" or \"mode2\", got \"$ic\"")
    end
    d0 = bc == :fixed ? diff([0.0; q0; 0.0]) : [q0[mod1(i+1, N)] - q0[i] for i in 1:N]
    strain0 = α * maximum(abs, d0)
    Q, V, t, _, _ = solve_fput(sp, q0, v0, (0.0, T), DT; saveat = range(0, T; length = NSAVE + 1))

    H   = [hamiltonian(Q[:, i], V[:, i], k, bc) for i in axes(Q, 2)]
    rel = abs.(H .- H[1]) ./ H[1]
    dec = length(rel) ÷ 10
    ok  = all(isfinite, Q) && maximum(rel) < 1e-5 && maximum(rel[end-dec+1:end]) < 3maximum(rel[2:dec])

    # Same projection as FPUTAnalysis.compute_modal_energies, on all columns at once
    # (column-by-column with reduce(hcat, …) over a generator is O(nt²) in copies).
    sm       = sqrt.(m)
    E_modal  = 0.5 .* ((Vm' * (sm .* V)) .^ 2 .+ freq .^ 2 .* (Vm' * (sm .* Q)) .^ 2)
    E_dimer  = dimer_energies(Q, V, 1 + dk, bc)

    # Optical share = modes above the gap. FBC has two in-gap modes at ω = √2
    # (N/2−1 acoustic, 2 gap, N/2−1 optical); they are counted in neither band.
    opt = bc == :fixed ? (N÷2+2:N) : (N÷2+1:N)
    opt_share = sum(E_modal[opt, end]) / sum(E_modal[:, end])
    # Share of the energy still in what was excited at t = 0.
    loc1 = ic == "dimer" ? (bc == :fixed ? E_dimer[2, end] + E_dimer[3, end] : E_dimer[1, end] + E_dimer[2, end]) /
                           sum(E_dimer[:, end]) :
                           E_modal[IC_MODE[ic], end] / sum(E_modal[:, end])

    (; dk, t_cyc = t .* ω_CLOCK ./ 2π, ok, dH = maximum(rel),
       S_modal = spectral_entropy(E_modal, SMOOTH_DELTA),
       S_dimer = spectral_entropy(E_dimer, SMOOTH_DELTA),
       opt_share, loc1, strain0,
       bw_opt = freq[end] - freq[first(opt)])
end

function main()
    N, E_in, cycles = Int(arg(1, 64)), arg(2, 0.445), arg(3, 1e5)
    deltas = length(ARGS) >= 4 ? parse.(Float64, split(ARGS[4], ",")) : DEFAULT_DELTAS
    ic     = length(ARGS) >= 5 ? ARGS[5] : "dimer"
    bc     = length(ARGS) >= 6 ? Symbol(ARGS[6]) : :periodic
    tag    = "N$(N)_dk" * join(deltas, "-") * (ic == "dimer" ? "" : "_$(ic)") * (bc == :fixed ? "_fbc" : "")
    T = cycles * 2π / ω_CLOCK
    @printf("%s, N=%d, α=%g, ic=%s with E=%g (ε=%.3g), T=%.3e (%.0e cycles of ω=2)\n\n",
            bc == :fixed ? "FBC" : "PBC", N, α, ic, E_in, E_in/N, T, cycles)
    @printf("%-6s %-5s %-9s %-9s %-10s %-10s %-10s %-10s %-10s %s\n",
            "Δκ", "ok", "|ΔH|/H", "α|d|(0)", "opt.bandw", "S_mod(0)", "S_mod(end)", "S_dim(end)", "E_opt/E",
            ic == "dimer" ? "E in dimer 1 (end)" : "E in $(ic) (end)")

    runs = map(deltas) do dk
        r = run(dk, N, E_in, T, ic, bc)
        @printf("%-6g %-5s %-9.2e %-9.3f %-10.3e %-10.3f %-10.3f %-10.3f %-10.3f %.3f\n",
                dk, r.ok, r.dH, r.strain0, r.bw_opt, r.S_modal[2], r.S_modal[end], r.S_dimer[end], r.opt_share, r.loc1)
        r
    end
    all(r -> r.ok, runs) || (println("\nA run failed the stability/energy check."); exit(1))
    @printf("\nreferences: ln N = %.3f, ln(N/2) = %.3f (whole optical band only)\n", log(N), log(N ÷ 2))

    jldsave(joinpath(mkpath(joinpath(@__DIR__, "data")), "dimer_approach_$(tag).jld2"); runs, N, E_in)

    save_fig(make_figure(runs, N, E_in, ic, cycles), joinpath(@__DIR__, "figures"),
             "dimer_approach_entropy_$(tag)"; exts = (".png",))
end

const IC_TITLE = Dict("dimer" => "Only dimer 1 stretched", "mode2" => "Mode-2 initial condition",
                      "mode1" => "Mode-1 initial condition")

"""Initial condition from a file name: `_mode1` / `_mode2` suffix, else the dimer."""
function ic_from_name(path)
    i = findfirst(ic -> occursin("_$(ic)", basename(path)), ("mode1", "mode2"))
    isnothing(i) ? "dimer" : ("mode1", "mode2")[i]
end

"""
    make_figure(runs, N, E_in, ic, cycles; normalize=false)

Two panels, one initial condition. The figure title states the initial condition;
each panel title names the set the energy is shared over, which is the only thing
that differs between them: the N normal modes, or the N dimer coordinates (N/2
relative + N/2 centre-of-mass, hence the ln N ceiling). Run parameters go in the
caption or the file name, not the figure.

`normalize=true` plots S̄/ln N, the fraction of the equipartition ceiling (1 ⇔ the
energy is shared evenly over all N elements), bounded in [0,1] and comparable
across N. The second reference line, ln(N/2)/ln N, is the ceiling when the energy
stays in one band.
"""
function make_figure(runs, N, E_in, ic, cycles; normalize::Bool = false)
    apply_style!()
    scale = normalize ? 1 / log(N) : 1.0
    panel(field, ttl) = begin
        p = logplot(; xlabel = L"t\ \mathrm{(cycles\ of\ }\omega=2)",
                    ylabel = normalize ? L"\bar{S}(t)\,/\,\ln N" : L"\bar{S}(t)",
                    decades = 0:Int(log10(cycles)), ylims = scale .* (-0.05log(N), 1.08log(N)),
                    legend = :bottomright, size = (700, 520), title = ttl, titlefontsize = 14)
        draw!(p, [Curve(logdownsample(r.t_cyc[2:end], scale .* getfield(r, field)[2:end])...,
                        latexstring("\\Delta\\kappa = $(r.dk)"); cyclic(j; lw = 2.2)...)
                  for (j, r) in enumerate(runs)])
        guide_hline!(p, scale * log(N);
                     annotation = normalize ? L"\mathit{equipartition}" : L"\ln N", side = :left)
        guide_hline!(p, scale * log(N ÷ 2);
                     annotation = normalize ? L"\ln(N/2)\,/\,\ln N" : L"\ln(N/2)", side = :left)
        p
    end
    plot(panel(:S_modal, "Normal modes"), panel(:S_dimer, "Dimer coordinates");
         layout = (1, 2), size = (1400, 580), plot_title = IC_TITLE[ic], plot_titlefontsize = 15,
         top_margin = 4Plots.mm)
end

"""
    replot(path; normalize=false)

Regenerate a figure from a saved `dimer_approach_*.jld2` without re-simulating.
The initial condition is read from the file name (`_mode2` suffix), the time span
from the saved clock.
"""
function replot(path; normalize::Bool = false)
    d  = load(path)
    ic = ic_from_name(path)
    cycles = exp10(round(log10(last(d["runs"][1].t_cyc))))
    name = replace(splitext(basename(path))[1], "dimer_approach_" => "dimer_approach_entropy_") *
           (normalize ? "_norm" : "")
    save_fig(make_figure(d["runs"], d["N"], d["E_in"], ic, cycles; normalize),
             joinpath(@__DIR__, "figures"), name; exts = (".png",))
end

"""
    merge_runs(paths...) -> String

Combine saved sweeps that share N, E and initial condition into one file, sorted by
Δκ, and return its path. A Δκ present in several inputs is kept once (first file wins),
after checking that the copies really are the same run — the simulation is
deterministic, so a mismatch means the inputs are not comparable.
"""
function merge_runs(paths...)
    ds = load.(paths)
    N, E_in = ds[1]["N"], ds[1]["E_in"]
    all(d -> d["N"] == N && d["E_in"] == E_in, ds) || error("Cannot merge: N or E differ between files")
    ics = unique(ic_from_name.(paths))
    fbc = unique(occursin("_fbc", basename(p)) for p in paths)
    length(ics) == 1 && length(fbc) == 1 || error("Cannot merge: mixed initial or boundary conditions")

    runs = Any[]
    for d in ds, r in d["runs"]
        i = findfirst(x -> isapprox(x.dk, r.dk), runs)
        if isnothing(i)
            push!(runs, r)
        else
            gap = max(maximum(abs, runs[i].S_dimer .- r.S_dimer), maximum(abs, runs[i].S_modal .- r.S_modal))
            gap < 1e-9 || error("Δκ=$(r.dk) appears twice with different results (max |ΔS̄| = $gap)")
            @printf("  Δκ=%-5g repeated, copies identical (max |ΔS̄| = %.1e): kept once\n", r.dk, gap)
        end
    end
    sort!(runs, by = r -> r.dk)

    tag = "N$(N)_dk" * join([r.dk for r in runs], "-") * (only(ics) == "dimer" ? "" : "_$(only(ics))") * (only(fbc) ? "_fbc" : "")
    out = joinpath(@__DIR__, "data", "dimer_approach_$(tag).jld2")
    jldsave(out; runs, N, E_in)
    println("merged Δκ = $([r.dk for r in runs]) → $out")
    out
end

if abspath(PROGRAM_FILE) == @__FILE__
    # `--replot [--norm] a.jld2 b.jld2 …` regenerates figures from saved runs.
    if !isempty(ARGS) && ARGS[1] == "--replot"
        rest = ARGS[2:end]
        normalize = !isempty(rest) && rest[1] == "--norm"
        foreach(f -> replot(f; normalize), normalize ? rest[2:end] : rest)
    elseif !isempty(ARGS) && ARGS[1] == "--merge"
        replot(merge_runs(ARGS[2:end]...); normalize = true)
    else
        main()
    end
end
