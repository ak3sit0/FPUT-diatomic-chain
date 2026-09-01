"""
    case_setup.jl

Turns `(spec, N, param, delta)` into a fully specified initial value problem.

This is the seam of the compute scripts: everything above it is config parsing,
everything below it is integration and reduction. It absorbs four blocks that
used to be copied verbatim across the five scripts — the `SystemParams` ternary,
the normal modes, the band selection, and the initial condition.

Two roles that used to share one variable are now separate fields:

- `excited`  — which mode(s) receive the energy.
- `ref_mode` — whose frequency defines `scaled_t` **and** the `TMAX` budget.

They coincide under `init_type="mode"`, which is why they were conflated. Under
`band_ensemble` there is no single excited mode, and `ref_mode` is the band's
lowest mode. Keeping them apart means nothing gets silently overwritten.
"""
module CaseSetup

using LinearAlgebra
using ..FPUTCore
using ..Experiment

export Case, build_case, initial_condition

"""
Fully specified case. `budget` carries the per-N `(E_total, TMAX, T_block, save_every)`.
"""
struct Case
    sp::SystemParams
    N::Int
    param::Float64
    delta::Float64
    freq::Vector{Float64}
    V::Matrix{Float64}
    m::Vector{Float64}
    excited::Vector{Int}        # who receives the energy (one entry under "mode")
    ref_mode::Int               # sets scaled_t and the TMAX budget
    omega_ref::Float64
    k_ac::Vector{Int}           # acoustic band, for interband diagnostics
    k_opt::Vector{Int}          # optical band
    budget::NamedTuple
end

"""
    system_params(spec, N, param, delta) -> SystemParams

`delta` is a spring or a mass disorder depending on `system_type`; `param` is the
cubic or quartic coupling depending on `nonlinear`. This dispatch-on-string used
to be repeated verbatim in four scripts.
"""
system_params(spec::SweepSpec, N::Integer, param::Real, delta::Real) = SystemParams(
    N,
    spec.system_type == "springs" ? Float64(delta) : 0.0,
    spec.system_type == "masses"  ? Float64(delta) : 0.0,
    spec.nonlinear == :alpha ? Float64(param) : 0.0,
    spec.nonlinear == :beta  ? Float64(param) : 0.0,
    spec.boundary,
)

"""
    bands(spec, N, freq) -> (excited, ref_mode, k_ac, k_opt)

Band indices and the two mode roles. Under `band_ensemble` the whole band is
excited and `ref_mode` is its lowest mode; an `init_mode` given there does not
apply and is reported rather than dropped in silence.
"""
function bands(spec::SweepSpec, N::Integer, freq::Vector{Float64})
    if spec.init_type == "band_ensemble"
        k_band = derive_band_indices(spec.branch, N, freq, spec.boundary;
                                     k_band_start = spec.k_band_start,
                                     k_band_end   = spec.k_band_end)
        other  = spec.branch == "acoustic" ? "optical" : "acoustic"
        k_opt  = spec.boundary == :fixed ? Int[] :
                 derive_band_indices(other, N, freq, spec.boundary)
        k_ac   = spec.branch == "acoustic" ? k_band : k_opt

        isnothing(spec.init_mode) || @info(
            "init_mode is ignored with init_type=\"band_ensemble\" (there is no single " *
            "excited mode); the reference frequency comes from the band's lowest mode",
            init_mode = spec.init_mode, ref_mode = first(k_band))

        return k_band, first(k_band), k_ac, isempty(k_opt) ? Int[] : k_opt
    else  # "mode"
        target = excited_mode(spec)
        return [target], target, collect(1:N÷2), collect(N÷2+1:N)
    end
end

"""
    build_case(spec, N, param, delta) -> Case

Assemble the case. `find_normal_modes` already returns ascending frequencies
(LAPACK `syevr` on a `Symmetric` matrix, and `sqrt` is monotone), so the
`sortperm(freq)` that used to follow every call was a no-op and is gone.
"""
function build_case(spec::SweepSpec, N::Integer, param::Real, delta::Real)
    sp      = system_params(spec, N, param, delta)
    k, m    = make_system(sp)
    freq, V = find_normal_modes(k, m, spec.boundary)

    excited, ref_mode, k_ac, k_opt = bands(spec, N, freq)
    ω_ref  = ref_frequency(freq, ref_mode)     # errors on the PBC ω≈0 translation
    budget = resolve_budget(spec, N, ω_ref)

    Case(sp, N, Float64(param), Float64(delta), freq, V, m,
         excited, ref_mode, ω_ref, k_ac, k_opt, budget)
end

"""
    initial_condition(case, spec, seed) -> (q0, v0)

Kept separate from `build_case` because an ensemble needs a different draw per
realization while the rest of the case is shared. Under `init_type="mode"` the
condition is deterministic and `seed` is ignored.
"""
function initial_condition(case::Case, spec::SweepSpec, seed::Integer)
    if spec.init_type == "band_ensemble"
        return band_phase_ic(case.excited, case.budget.E_total, case.N,
                             case.freq, case.V, case.m, seed)
    else
        target    = only(case.excited)
        U         = Diagonal(1.0 ./ sqrt.(case.m)) * case.V
        amplitude = sqrt(2 * case.budget.E_total) / case.freq[target]
        return amplitude .* U[:, target], zeros(case.N)
    end
end

end # module CaseSetup
