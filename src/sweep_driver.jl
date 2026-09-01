"""
    sweep_driver.jl

The `main` loop that every compute script repeated: run a list of tasks, survive
individual failures, keep what worked, save it with the config path.

Deliberately dumb about physics. It does not decide the task order (cost
heuristics belong to the caller, which knows that cost grows as N·TMAX ∝ N²) nor
the parallelism policy — that is an argument, because the right answer differs:
trajectories parallelize over cases, an ensemble parallelizes over realizations
*inside* each case and so must run its cases serially.
"""
module SweepDriver

using JLD2, Dates, Base.Threads

export run_sweep

"""
    run_sweep(f, tasks; outfile, config_path, parallel=:threads, label="case") -> String

Apply `f` to each task, save the non-`nothing` results to `outfile`, return its
path. A task that throws is reported with its backtrace and dropped, so one bad
parameter combination cannot lose a whole sweep.

`parallel`:
- `:threads` — `@threads` over tasks (uniform cost).
- `:dynamic` — `@threads :dynamic`, for uneven cost; pre-sort `tasks` yourself.
- `:serial`  — one task at a time, for when `f` is itself threaded.
"""
function run_sweep(f, tasks; outfile::AbstractString, config_path::AbstractString,
                   parallel::Symbol = :threads, label::AbstractString = "case")
    n = length(tasks)
    n > 0 || error("Empty task list")
    results = Vector{Any}(undef, n)

    attempt(i) = try
        f(tasks[i])
    catch e
        println("[$label $i] ERROR: $e\n$(sprint(showerror, e, catch_backtrace()))")
        nothing
    end

    if parallel === :threads
        @threads for i in 1:n
            results[i] = attempt(i)
        end
    elseif parallel === :dynamic
        @threads :dynamic for i in 1:n
            results[i] = attempt(i)
        end
    elseif parallel === :serial
        for i in 1:n
            results[i] = attempt(i)
        end
    else
        error("parallel must be :threads, :dynamic or :serial; got :$parallel")
    end

    valid = filter(!isnothing, results)
    isempty(valid) && @warn "No task produced a result; saving an empty file" outfile

    mkpath(dirname(outfile))
    jldsave(outfile; results = valid, config = config_path)
    println("Saved: $outfile  ($(length(valid))/$n valid cases)")
    outfile
end

"""
    sweep_tasks(spec) -> Vector{@NamedTuple{N::Int, param::Float64, delta::Float64}}

The full `N × param × Δ` grid, heaviest first so a dynamic schedule balances:
cost grows as `N·TMAX`, and `TMAX ∝ N` whenever `scaled_t_max` is set.
"""
function sweep_tasks(spec)
    tasks = [(; N = N, param = p, delta = d)
             for N in spec.N_values, p in spec.param_values, d in spec.delta_values]
    sort!(vec(tasks), by = t -> -t.N)
end

export sweep_tasks

"""Default output path: `<base_dir>/<stem>_<today>.jld2`."""
output_path(spec, stem::AbstractString) =
    joinpath(spec.base_dir, "$(stem)_$(Dates.today()).jld2")

export output_path

end # module SweepDriver
