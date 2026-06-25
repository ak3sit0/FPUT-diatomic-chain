"""
    check_hpc_status.jl

Revisa el estado de los jobs HPC multi-N y genera un reporte.
Verifica logs, archivos .SUCCESS/.FAILED, y si existen outputs JLD2.

Usage:
  julia --project=. scripts/check_hpc_status.jl

Output:
  - Tabla de estado de cada N
  - Detalles de errores si los hay
  - Comandos para ver logs completos
"""

using Printf

function check_job_status(N::Int)
    """Verifica el estado de un job N específico."""
    # Buscar archivos
    logs = readdir("results/logs", join=true) |>
           x -> filter(f -> match(Regex("ensemble_N$(N).*\\.log"), f) !== nothing, x)

    success = readdir("results/logs", join=true) |>
              x -> filter(f -> match(Regex("ensemble_N$(N).*\\.SUCCESS"), f) !== nothing, x)

    failed = readdir("results/logs", join=true) |>
             x -> filter(f -> match(Regex("ensemble_N$(N).*\\.FAILED"), f) !== nothing, x)

    # Buscar output JLD2
    output_dir = "results/data/ensemble_N_sweep_N$(N)"
    output_exists = isdir(output_dir)
    output_size = 0
    if output_exists
        jld2_files = filter(f -> endswith(f, ".jld2"), readdir(output_dir))
        if !isempty(jld2_files)
            output_size = sum(filesize(joinpath(output_dir, f)) for f in jld2_files)
        end
    end

    # Determinar estado
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
            output_exists=output_exists, output_size=output_size)
end

function main()
    N_values = [32, 64, 128, 256]

    println()
    println("=" ^ 80)
    println("HPC JOB STATUS REPORT")
    println("=" ^ 80)
    println()

    # Tabla de resumen
    println(@sprintf "%-8s | %-20s | %-10s | %-15s", "N", "Status", "Output", "Size")
    println("-" ^ 60)

    results = Dict()
    for N in N_values
        r = check_job_status(N)
        results[N] = r

        output_str = r.output_exists ? "✓ Exists" : "✗ Missing"
        size_str = r.output_exists ? "$(round(r.output_size / 1e9; digits=2)) GB" : "—"

        println(@sprintf "%-8d | %-20s | %-10s | %-15s",
                N, r.status, output_str, size_str)
    end

    println()
    println("=" ^ 80)
    println()

    # Detalles de cada job
    for N in N_values
        r = results[N]

        if r.status == "✓ SUCCESS"
            println("[$N] ✓ Completed successfully")
            if !isempty(r.logs)
                println("     Log: $(basename(r.logs[1]))")
            end
            if r.output_exists
                println("     Output: results/data/ensemble_N_sweep_N$(N)/")
            end
        elseif r.status == "✗ FAILED"
            println("[$N] ✗ FAILED - Ver detalles:")
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
            println("[$N] ⏳ Still running or in queue")
            if !isempty(r.logs)
                println("     Log: $(r.logs[1])")
                println("     Ver progreso en tiempo real:")
                println("       tail -f $(r.logs[1])")
            end
        else
            println("[$N] ❓ No logs found yet")
            println("     Comando para enviar:")
            println("       qsub jobs/ensemble_N$(N).pbs")
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
        for N in N_values
            if results[N].status == "✗ FAILED"
                println("  qsub jobs/ensemble_N$(N).pbs")
            end
        end
    elseif running_count > 0
        println("ℹ️  $(running_count) jobs still running. Check again later.")
        println()
        println("To monitor in real-time:")
        println("  qstat")
        println("  tail -f results/logs/ensemble_N*.log")
    elseif success_count == length(N_values)
        println("✅ ALL JOBS COMPLETED SUCCESSFULLY!")
        println()
        println("Next steps:")
        println("  1. Analyze results:")
        println("     julia --project=. examples/plot_thermalization_time.jl \\")
        println("       results/data/ensemble_N_sweep_N32/ensemble_results_*.jld2")
        println()
        println("  2. Check output sizes (should be >100MB each):")
        println("     du -h results/data/ensemble_N_sweep_N*/")
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
