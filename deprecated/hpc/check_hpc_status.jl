"""
    check_hpc_status.jl

Revisa el estado de los jobs HPC y genera un reporte.
Verifica logs, archivos .SUCCESS/.FAILED, y si existen outputs JLD2.

Descubre los jobs a partir de los .pbs generados en jobs/ (por generate_hpc_jobs.jl
o generate_pbs_nsweep.jl indistintamente) en vez de asumir un naming/N_values fijo,
así funciona con cualquiera de los dos generadores.

Usage:
  julia --project=. scripts/hpc/check_hpc_status.jl
"""

using Printf
using TOML

"Nombres de job (basename sin extensión) de todos los .pbs en jobs/."
function discover_jobs(jobs_dir::String="jobs")
    isdir(jobs_dir) || return String[]
    pbs_files = filter(f -> endswith(f, ".pbs"), readdir(jobs_dir))
    sort([splitext(f)[1] for f in pbs_files])
end

"Lee output.base_dir del .toml homónimo del job, si existe."
function output_dir_for_job(jobname::String, jobs_dir::String="jobs")
    toml_path = joinpath(jobs_dir, "$(jobname).toml")
    isfile(toml_path) || return nothing
    cfg = TOML.parsefile(toml_path)
    get(get(cfg, "output", Dict{String,Any}()), "base_dir", nothing)
end

"""
    check_job_status(jobname; jobs_dir="jobs", logs_dir="results/logs")

Verifica el estado de un job específico buscando `<jobname>_<PBS_JOBID>.{log,SUCCESS,FAILED}`
en `logs_dir`, y el directorio de salida declarado en `jobs/<jobname>.toml`.
"""
function check_job_status(jobname::String; jobs_dir::String="jobs", logs_dir::String="results/logs")
    esc = replace(jobname, r"([.^$|()\[\]{}*+?\\])" => s"\\\1")
    pattern = Regex("^$(esc)_.*\\.(log|SUCCESS|FAILED)\$")
    files = isdir(logs_dir) ? filter(f -> match(pattern, f) !== nothing, readdir(logs_dir)) : String[]

    logs    = joinpath.(logs_dir, filter(f -> endswith(f, ".log"), files))
    success = joinpath.(logs_dir, filter(f -> endswith(f, ".SUCCESS"), files))
    failed  = joinpath.(logs_dir, filter(f -> endswith(f, ".FAILED"), files))

    output_dir = output_dir_for_job(jobname, jobs_dir)
    output_exists = !isnothing(output_dir) && isdir(output_dir)
    output_size = 0
    if output_exists
        jld2_files = filter(f -> endswith(f, ".jld2"), readdir(output_dir))
        if !isempty(jld2_files)
            output_size = sum(filesize(joinpath(output_dir, f)) for f in jld2_files)
        end
    end

    status = if !isempty(success)
        "✓ SUCCESS"
    elseif !isempty(failed)
        "✗ FAILED"
    elseif !isempty(logs)
        "⏳ RUNNING/QUEUED"
    else
        "❓ NO LOGS"
    end

    return (status=status, logs=logs, success=success, failed=failed,
            output_dir=output_dir, output_exists=output_exists, output_size=output_size)
end

function main()
    jobnames = discover_jobs()
    if isempty(jobnames)
        println("No se encontraron .pbs en jobs/. Generar jobs primero con")
        println("  julia --project=. scripts/hpc/generate_hpc_jobs.jl <config>")
        println("  julia --project=. scripts/hpc/generate_pbs_nsweep.jl <config>")
        return
    end

    println()
    println("=" ^ 80)
    println("HPC JOB STATUS REPORT")
    println("=" ^ 80)
    println()

    println(@sprintf("%-28s | %-20s | %-10s | %-15s", "Job", "Status", "Output", "Size"))
    println("-" ^ 80)

    results = Dict{String,Any}()
    for jobname in jobnames
        r = check_job_status(jobname)
        results[jobname] = r

        output_str = r.output_exists ? "✓ Exists" : "✗ Missing"
        size_str = r.output_exists ? "$(round(r.output_size / 1e9; digits=2)) GB" : "—"

        println(@sprintf("%-28s | %-20s | %-10s | %-15s",
                jobname, r.status, output_str, size_str))
    end

    println()
    println("=" ^ 80)
    println()

    # Detalles de cada job
    for jobname in jobnames
        r = results[jobname]

        if r.status == "✓ SUCCESS"
            println("[$jobname] ✓ Completed successfully")
            if !isempty(r.logs)
                println("     Log: $(basename(r.logs[1]))")
            end
            if r.output_exists
                println("     Output: $(r.output_dir)/")
            end
        elseif r.status == "✗ FAILED"
            println("[$jobname] ✗ FAILED - Ver detalles:")
            if !isempty(r.failed)
                failed_file = r.failed[1]
                println("     Failed marker: $(basename(failed_file))")
                if !isempty(r.logs)
                    log_file = r.logs[1]
                    println("     Log: $(log_file)")
                    println("\n     --- LAST 30 LINES OF LOG ---")
                    try
                        lines = readlines(log_file)
                        for line in lines[max(1, length(lines)-30):end]
                            println("     $line")
                        end
                    catch e
                        println("     Error reading log: $e")
                    end
                    println("     --- END OF LOG ---\n")
                    println("     Para ver el log completo:")
                    println("       cat $(log_file)")
                end
            end
        elseif r.status == "⏳ RUNNING/QUEUED"
            println("[$jobname] ⏳ Still running or in queue")
            if !isempty(r.logs)
                println("     Log: $(r.logs[1])")
                println("     Ver progreso en tiempo real:")
                println("       tail -f $(r.logs[1])")
            end
        else
            println("[$jobname] ❓ No logs found yet")
            println("     Comando para enviar:")
            println("       qsub jobs/$(jobname).pbs")
        end

        println()
    end

    # Sumario final
    success_count = sum(1 for r in values(results) if r.status == "✓ SUCCESS")
    failed_count = sum(1 for r in values(results) if r.status == "✗ FAILED")
    running_count = sum(1 for r in values(results) if r.status == "⏳ RUNNING/QUEUED")
    pending_count = sum(1 for r in values(results) if r.status == "❓ NO LOGS")

    println("=" ^ 80)
    println("SUMMARY")
    println("=" ^ 80)
    println("  ✓ Success:  $success_count")
    println("  ✗ Failed:   $failed_count")
    println("  ⏳ Running:  $running_count")
    println("  ❓ Pending:  $pending_count")
    println()

    if failed_count > 0
        println("⚠️  ALERT: $failed_count jobs failed. Check logs above.")
        println()
        println("To resubmit a failed job:")
        for jobname in jobnames
            if results[jobname].status == "✗ FAILED"
                println("  qsub jobs/$(jobname).pbs")
            end
        end
    elseif running_count > 0
        println("ℹ️  $(running_count) jobs still running. Check again later.")
        println()
        println("To monitor in real-time:")
        println("  qstat")
        println("  tail -f results/logs/*.log")
    elseif success_count == length(jobnames)
        println("✅ ALL JOBS COMPLETED SUCCESSFULLY!")
        println()
        println("Next steps:")
        println("  1. Analyze results (see each job's output_dir above)")
        println("  2. Check output sizes (should be >100MB each):")
        println("     du -h results/data/*/")
    end

    println()
    println("=" ^ 80)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isdir("results/logs")
        main()
    else
        println("Error: results/logs/ not found. Run jobs first.")
    end
end
