# Code for generating animations of the FPUT-System
# Usage: julia --project=. examples/make_animation.jl

include("../src/parameters.jl"); using Main.Parameters
include("../src/dynamical_matrix.jl"); using Main.DynamicalMatrix
include("../src/fput_equations.jl"); using Main.FPUTEquations
include("../src/integrator.jl"); using Main.Integrator
include("../src/energy_analysis.jl"); using Main.EnergyAnalysis
include("../src/plotting.jl"); using Main.Plotting

using LinearAlgebra, Plots

# ============================================================================
# CONFIGURATION
# ============================================================================

function build_animation_config()
    return (
        N           = 32,           # chain length — change as needed
        init_mode   = 1,
        energy      = 0.45,         # initial energy per mode
        TMAX        = 10000.0,
        save_every  = 200,
        DeltaK      = 0.0,
        DeltaM      = 0.0,
        boundary    = :periodic,    # :periodic | :fixed
        start_time  = 0.0,
        end_time    = 1000.0,
        fps         = 20,
        marker_size = 6,
        color_odd   = RGB(0.12, 0.47, 0.71),    # Blue
        color_even  = RGB(1.0, 0.5, 0.05),      # Orange
        color_mono  = RGB(0.12, 0.47, 0.71),    # Blue
        line_color  = :gray,
        line_width  = 1.0,
    )
end

# ============================================================================
# SYSTEM BUILDER
# ============================================================================

struct FPUTSystem
    k::Vector{Float64}
    m::Vector{Float64}
    frequencies::Vector{Float64}
    U::Matrix{Float64}
    is_diatomic::Bool
    N::Int
end

function build_fput_system(N, DeltaK, DeltaM, boundary)
    k, m = Parameters.make_k_m(N, DeltaK, DeltaM, boundary)
    frequencies, normal_matrix = DynamicalMatrix.find_normal_modes(k, m)
    U = Diagonal(sqrt.(1 ./ m)) * normal_matrix
    is_diatomic = (length(unique(m)) > 1) || (length(unique(k)) > 1)
    FPUTSystem(k, m, frequencies, U, is_diatomic, N)
end

# ============================================================================
# INITIAL CONDITIONS
# ============================================================================

function make_initial_conditions(sys::FPUTSystem, init_mode, energy)
    amplitude = sqrt(2 * energy) / sys.frequencies[init_mode]
    q0 = amplitude * sys.U[:, init_mode]
    v0 = zeros(sys.N)
    (q0, v0)
end

# ============================================================================
# SIMULATION RUNNER
# ============================================================================

struct SimulationResult
    q::Matrix{Float64}
    v::Matrix{Float64}
    times::Vector{Float64}
    ylim::Tuple{Float64, Float64}
end

function run_fput_simulation(sys, ic, cfg)
    q0, v0 = ic
    alpha = 0.0
    beta = 0.0
    forces = zeros(sys.N)
    
    tspan = (cfg.start_time, cfg.end_time)
    
    sol = Integrator.run_simulation(v0, q0, tspan, sys.k, sys.m, alpha, beta, 
                                    Parameters.DT, forces; save_every=cfg.save_every)
    
    Nt = length(sol.u)
    q = zeros(sys.N, Nt)
    v = zeros(sys.N, Nt)
    for (idx, u) in enumerate(sol.u)
        v[:, idx] = u.x[1]
        q[:, idx] = u.x[2]
    end
    
    times = copy(sol.t)
    max_q = maximum(abs.(q))
    ylim = (-1.2 * max_q, 1.2 * max_q)
    
    SimulationResult(q, v, times, ylim)
end

# ============================================================================
# FRAME RENDERING
# ============================================================================

function make_renderer_diatomic(odd_idx, even_idx, cfg)
    return (q_frame, i) -> begin
        plt = plot(1:cfg.N, q_frame,
                   linecolor=cfg.line_color, linewidth=cfg.line_width,
                   legend=false, grid=false, axis=false, framestyle=:none,
                   markerstrokewidth=0)
        scatter!(plt, odd_idx, q_frame[odd_idx],
                 color=cfg.color_odd, markersize=cfg.marker_size,
                 markerstrokecolor=:transparent, markerstrokewidth=0)
        scatter!(plt, even_idx, q_frame[even_idx],
                 color=cfg.color_even, markersize=cfg.marker_size,
                 markerstrokecolor=:transparent, markerstrokewidth=0)
        plt
    end
end

function make_renderer_monatomic(cfg)
    return (q_frame, i) -> begin
        plot(1:cfg.N, q_frame,
             linecolor=cfg.line_color, linewidth=cfg.line_width,
             legend=false, grid=false, axis=false, framestyle=:none,
             markerstrokewidth=0)
        scatter!(1:cfg.N, q_frame, color=cfg.color_mono, 
                 markersize=cfg.marker_size, markerstrokecolor=:transparent,
                 markerstrokewidth=0)
    end
end

# ============================================================================
# ANIMATION BUILDER
# ============================================================================

function build_animation(result::SimulationResult, sys::FPUTSystem, cfg)
    Nt = length(result.times)
    println("Generating animation with $Nt frames...")
    
    if sys.is_diatomic
        odd_idx = 1:2:cfg.N
        even_idx = 2:2:cfg.N
        renderer_fn = make_renderer_diatomic(odd_idx, even_idx, cfg)
    else
        renderer_fn = make_renderer_monatomic(cfg)
    end
    
    anim = @animate for i in 1:Nt
        q_frame = result.q[:, i]
        plt = renderer_fn(q_frame, i)
        ylims!(plt, result.ylim)
        plt
    end
    
    return anim
end

# ============================================================================
# SAVING
# ============================================================================

function save_animation(anim, result::SimulationResult, cfg)
    output_dir = joinpath("results", "figures", "animations")
    mkpath(output_dir)
    
    timestr = "_t$(round(result.times[1], digits=2))_to_$(round(result.times[end], digits=2))"
    output_path = joinpath(output_dir, 
        "fput_animation_N$(cfg.N)_mode$(cfg.init_mode)_E$(round(cfg.energy, digits=3))$timestr.gif")
    
    println("Saving animation to $output_path...")
    gif(anim, output_path, fps=cfg.fps)
    println("Done!")
end

# ============================================================================
# MAIN ENTRY
# ============================================================================

function main()
    cfg = build_animation_config()
    println("Starting simulation for animation (N=$(cfg.N), mode=$(cfg.init_mode), E=$(cfg.energy), TMAX=$(cfg.TMAX))...")
    
    sys = build_fput_system(cfg.N, cfg.DeltaK, cfg.DeltaM, cfg.boundary)
    ic = make_initial_conditions(sys, cfg.init_mode, cfg.energy)
    result = run_fput_simulation(sys, ic, cfg)
    anim = build_animation(result, sys, cfg)
    save_animation(anim, result, cfg)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
