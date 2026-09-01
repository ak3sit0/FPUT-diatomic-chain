"""
    experiment.jl

Reads a sweep TOML once, into one `SweepSpec`, with **one name per concept**.

Two rules make this module worth existing:

1. It stores the *policy*, not the resolved numbers. `E_total`, `TMAX` and
   `save_every` all depend on N, and N varies inside a sweep — resolving them at
   parse time silently breaks any N sweep. `resolve_budget(spec, N, ω_ref)` does
   the per-N resolution.
2. The excited mode has a physical default per boundary, in one place with its
   reason: under PBC mode 1 is the k=0 uniform translation (ω≈0), so periodic
   starts at 2; fixed has no zero mode, so it starts at 1.
"""
module Experiment

using TOML

export SweepSpec, parse_spec, default_excited_mode, excited_mode, resolve_budget

"""
Parsed sweep configuration. `nothing` in a policy field means "use the other one":
`energy_density` xor `E_total`, `scaled_t_max` xor `TMAX`, `n_samples` xor `save_every`.
"""
struct SweepSpec
    # sweep axes
    N_values::Vector{Int}
    param_values::Vector{Float64}
    delta_values::Vector{Float64}
    # system
    boundary::Symbol
    system_type::String
    nonlinear::Symbol
    # initial condition
    init_type::String                     # "mode" | "band_ensemble"
    init_mode::Union{Int,Nothing}         # nothing → default_excited_mode(boundary)
    branch::String
    k_band_start::Union{Int,Nothing}
    k_band_end::Union{Int,Nothing}
    n_real::Int
    seed_base::Int
    # energy policy
    energy_density::Union{Float64,Nothing}
    E_total::Float64
    # time budget policy
    scaled_t_max::Union{Float64,Nothing}
    TMAX::Float64
    T_block::Float64
    n_blocks::Int
    DT::Float64
    # sampling policy
    n_samples::Union{Int,Nothing}
    save_every::Int
    downsample::Int
    # misc
    entropy_delta::Float64
    debug::Bool
    base_dir::String
end

_maybe(d, key, T) = haskey(d, key) ? T(d[key]) : nothing

"""
    parse_spec(path) -> SweepSpec

Parse a sweep TOML. Warns on `initial_condition`, the legacy spelling of
`init_type` — it was read into a variable that no script ever used, which is how
`init_mode` came to be ignored by `compute_trajectories.jl`.
"""
function parse_spec(path::AbstractString)
    isfile(path) || error("Config not found: $path")
    d    = TOML.parsefile(path)
    phys = d["physics"]
    sim  = d["simulation"]
    out  = d["output"]

    haskey(phys, "initial_condition") && !haskey(phys, "init_type") &&
        @warn "TOML uses 'initial_condition' (deprecated); interpreted as init_type" path

    init_type = get(phys, "init_type", get(phys, "initial_condition", "mode"))

    N_values = haskey(phys, "N_values") ? Vector{Int}(phys["N_values"]) : [Int(phys["N"])]
    isempty(N_values) && error("N_values empty in $path")

    SweepSpec(
        N_values,
        Float64.(phys["param_values"]),
        Float64.(phys["delta_values"]),
        Symbol(phys["boundary"]),
        String(phys["system_type"]),
        Symbol(phys["nonlinear"]),
        String(init_type),
        _maybe(phys, "init_mode", Int),
        String(get(phys, "branch", "acoustic")),
        _maybe(phys, "k_band_start", Int),
        _maybe(phys, "k_band_end", Int),
        Int(get(phys, "n_real", 1)),
        Int(get(phys, "seed_base", 42)),
        _maybe(phys, "energy_density", Float64),
        Float64(get(phys, "initial_energy", 0.445)),
        _maybe(sim, "scaled_t_max", Float64),
        Float64(get(sim, "TMAX", 0.0)),
        Float64(get(sim, "T_block", 0.0)),
        Int(get(sim, "n_blocks", 20)),
        Float64(get(sim, "DT", 0.05)),
        _maybe(sim, "n_samples", Int),
        Int(get(sim, "save_every", 1000)),
        Int(get(sim, "downsample", 1)),
        Float64(get(sim, "entropy_delta", 0.6)),
        Bool(get(sim, "debug", false)),
        String(out["base_dir"]),
    )
end

"""
    default_excited_mode(boundary) -> Int

Lowest physical mode. Periodic mode 1 is the uniform k=0 translation (ω≈0, does
not oscillate), so periodic starts at 2; fixed ends are clamped, so there is no
zero mode and mode 1 is already the lowest acoustic one.
"""
default_excited_mode(boundary::Symbol) = boundary == :fixed ? 1 : 2

"""
    excited_mode(spec) -> Int

Mode that receives the energy under `init_type="mode"`. Reports an explicit
`init_mode` that departs from the physical default — informational, not a
warning: departing is legitimate, being unable to see it was the problem.
"""
function excited_mode(spec::SweepSpec)
    default = default_excited_mode(spec.boundary)
    isnothing(spec.init_mode) && return default
    spec.init_mode != default &&
        @info "init_mode=$(spec.init_mode) explicit (the default for boundary=$(spec.boundary) would be $default)"
    spec.init_mode
end

"""
    resolve_budget(spec, N, ω_ref) -> (; E_total, TMAX, T_block, save_every)

Resolve the per-N quantities. Kept out of `parse_spec` on purpose:

- `energy_density` fixes ε ⇒ `E_total = ε·N`. Fixing `E_total` instead would make
  ε fall as 1/N, so large systems get progressively more linear.
- `scaled_t_max` fixes the observation window in cycles ⇒ `TMAX = scaled_t_max·2π/ω_ref`.
  Since `ω_ref ≈ 6.25/N`, that means `TMAX ∝ N`.
- `n_samples` fixes the number of stored samples ⇒ `modal_E ∝ N`, not `∝ N²`.
"""
function resolve_budget(spec::SweepSpec, N::Integer, ω_ref::Real)
    ω_ref > 0 || error("ω_ref=$ω_ref must be positive (use FPUTCore.ref_frequency)")

    E_total = isnothing(spec.energy_density) ? spec.E_total : spec.energy_density * N

    TMAX = isnothing(spec.scaled_t_max) ? spec.TMAX : spec.scaled_t_max * 2π / ω_ref
    TMAX > 0 || error("Config must give 'TMAX' or 'scaled_t_max' in [simulation]")

    save_every = isnothing(spec.n_samples) ? spec.save_every :
                 max(1, round(Int, TMAX / (spec.DT * spec.n_samples)))

    T_block = spec.T_block > 0 ? min(spec.T_block, TMAX) : TMAX / spec.n_blocks

    (; E_total, TMAX, T_block, save_every)
end

end # module Experiment
