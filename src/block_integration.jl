module BlockIntegration

using ..FPUTFastRunner
using ..FPUTAnalysis

export integrate_in_blocks

"""
    integrate_in_blocks(sp, q0, v0, N, freq, V, m, target_mode_idx;
                         TMAX, T_block, DT, save_every, downsample,
                         track_abs_time=false, debug=false, label="",
                         on_block=nothing, check_instability=nothing)
      -> (scaled_t, t_abs, modal_E, q_final, v_final, aborted)

Block-time-stepping integration loop shared by `compute_trajectories*.jl` and
`compute_ensemble*.jl`: integrates `sp` from `(q0,v0)` in chunks of `T_block`
up to `TMAX`, normalizes Q/V shape to N×nt, downsamples avoiding a duplicated
sample at block boundaries, and accumulates modal energies into one matrix.
Time is expressed in cycles of `freq[target_mode_idx]` (scaled_t = t·ω/2π).

- `on_block(modal_Eb, idx_ds)`: called once per non-empty block (after
  downsampling) for extra per-block accumulation (e.g. per-band energy sums).
- `check_instability(modal_Eb)`: called once per block; if it returns `true`,
  integration stops early and `aborted=true` is returned.
- `track_abs_time`: also accumulate the untransformed (absolute) time samples
  in `t_abs`; empty if `false`.
"""
function integrate_in_blocks(sp, q0, v0, N::Int, freq::Vector{Float64},
                              V::Matrix{Float64}, m::Vector{Float64},
                              target_mode_idx::Int;
                              TMAX::Float64, T_block::Float64, DT::Float64,
                              save_every::Int, downsample::Int,
                              track_abs_time::Bool=false, debug::Bool=false,
                              label::String="",
                              on_block=nothing, check_instability=nothing)
    q_cur = copy(q0)
    v_cur = copy(v0)
    T_total        = Float64[]
    T_abs_total    = Float64[]
    modal_E_blocks = Vector{Matrix{Float64}}()
    t_cur          = 0.0
    aborted        = false

    while t_cur < TMAX
        t_next = min(t_cur + T_block, TMAX)
        saveat = t_cur:save_every*DT:t_next

        Qb, Vb, Tb, _, _ = FPUTFastRunner.solve_fput(sp, q_cur, v_cur,
                                                       (t_cur, t_next), DT;
                                                       saveat=saveat)

        Tb_abs = track_abs_time ? Float64.(Tb) : Float64[]
        Tb = Tb .* freq[target_mode_idx] ./ (2π)

        if size(Qb,1) == N && size(Qb,2) >= 1
            Qmat = Float64.(Qb)
            Vmat = Float64.(Vb)
        elseif size(Qb,2) == N && size(Qb,1) >= 1
            Qmat = Float64.(Qb')
            Vmat = Float64.(Vb')
        else
            error("Unexpected Q/V shape: Q=$(size(Qb)), V=$(size(Vb))")
        end

        modal_Eb = FPUTAnalysis.compute_modal_energies(Qmat, Vmat, freq, V, m)

        if check_instability !== nothing && check_instability(modal_Eb)
            aborted = true
            break
        end

        if !isempty(T_total) && !isempty(Tb) && isapprox(T_total[end], Tb[1]; atol=1e-12, rtol=0)
            if length(Tb) > 1
                idx_ds = 2:downsample:length(Tb)
                append!(T_total, Float64.(Tb[idx_ds]))
                track_abs_time && append!(T_abs_total, Tb_abs[idx_ds])
                blk = Float64.(modal_Eb[:, idx_ds])
                size(blk,2) > 0 && push!(modal_E_blocks, blk)
                on_block !== nothing && on_block(modal_Eb, idx_ds)
            end
        else
            idx_ds = 1:downsample:length(Tb)
            append!(T_total, Float64.(Tb[idx_ds]))
            track_abs_time && append!(T_abs_total, Tb_abs[idx_ds])
            blk = Float64.(modal_Eb[:, idx_ds])
            size(blk,2) > 0 && push!(modal_E_blocks, blk)
            on_block !== nothing && on_block(modal_Eb, idx_ds)
        end

        q_cur .= Qmat[:, end]
        v_cur .= Vmat[:, end]
        t_cur  = t_next
        debug && println("[$label] t=$(round(t_cur; digits=3)) / $(TMAX)")
    end

    modal_E = isempty(modal_E_blocks) ? zeros(N, 0) : reduce(hcat, modal_E_blocks)

    return (scaled_t=T_total, t_abs=T_abs_total, modal_E=modal_E,
            q_final=q_cur, v_final=v_cur, aborted=aborted)
end

end # module BlockIntegration
