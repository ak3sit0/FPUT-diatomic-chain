module FPUTEquations

export fput_generalized

function fput_generalized(dv, v, q, p, t)
    k, m, alpha, beta, forces = p
    N = length(m)
    forces .= 0.0 # Reset the forces array

    # Compute forces with nonlinear terms and different boundary conditions

    if length(k) == N + 1 # Fixed boundary conditions
        for i in 2:N-1
            dl = q[i] - q[i-1]
            dr = q[i+1] - q[i]
            forces[i] = k[i+1]*dr - k[i]*dl + alpha*(k[i+1]*dr^2 - k[i]*dl^2) + beta*(k[i+1]*dr^3 - k[i]*dl^3)
        end
        forces[1] = k[2]*(q[2]-q[1]) - k[1]*q[1] + alpha*(k[2]*(q[2]-q[1])^2 - k[1]*q[1]^2) + beta*(k[2]*(q[2]-q[1])^3 - k[1]*q[1]^3)
        
        dl_N = q[N] - q[N-1]
        dr_N = -q[N]
        forces[N] = k[N+1]*dr_N - k[N]*dl_N + alpha*(k[N+1]*dr_N^2 - k[N]*dl_N^2) + beta*(k[N+1]*dr_N^3 - k[N]*dl_N^3)
    elseif length(k) == N # Periodic boundary conditions
        for i in 1:N
            #dl = q[i] - q[mod1(i-1, N)]
            #dr = q[mod1(i+1, N)] - q[i]
            #forces[i] = k[mod1(i+1, N)]*dr - k[i]*dl + alpha*(k[mod1(i+1, N)]*dr^2 - k[i]*dl^2) + beta*(k[mod1(i+1, N)]*dr^3 - k[i]*dl^3)

            # Índices con wrapping (periódicos)
            idx_prev = (i == 1) ? N : i - 1
            idx_next = (i == N) ? 1 : i + 1
            
            # Deformaciones relativas
            r_right = q[idx_next] - q[i]
            r_left  = q[i] - q[idx_prev]
            
            # Fuerza neta: F(derecha) - F(izquierda)
            forces[i] = (k[i]*r_right + alpha*r_right^2 + beta*r_right^3) - 
                   (k[idx_prev]*r_left + alpha*r_left^2 + beta*r_left^3)
        end
    else
        error("Inconsistent k and m lengths for boundary conditions.")
    end

    #for i in 2:N-1
    #    dl = q[i] - q[i-1]
    #    dr = q[i+1] - q[i]
    #    forces[i] = k[i+1]*dr - k[i]*dl + alpha*(k[i+1]*dr^2 - k[i]*dl^2) + beta*(k[i+1]*dr^3 - k[i]*dl^3)
    #end
    #forces[1] = -k[1]*q[1] + k[2]*(q[2] - q[1]) + alpha*(k[2]*(q[2]-q[1])^2 - k[1]*q[1]^2) + beta*(k[2]*(q[2]-q[1])^3 - k[1]*q[1]^3)
    #forces[N] = -k[N+1]*q[N] + k[N]*(q[N-1] - q[N]) + alpha*(-k[N]*(q[N-1]-q[N])^2 + k[N+1]*q[N]^2) + beta*(-k[N]*(q[N-1]-q[N])^3 - k[N+1]*q[N]^3)
    
    dv .= forces ./ m
end

end