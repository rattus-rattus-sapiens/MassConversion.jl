using MassConversion

function main(reps::Int64, blocksize::Int64)
    println("Initialising")
    t_span = 0.0:0.01:10
    D_init = [0]
    C_init = [1000]
    λ_reac = [1e0, 0.0]
    λ_tran = 1.0
    R_mat = [
        -1 0
        0 0
        ]
        θ = [(150, 200)]
        
        @inline function F!(dxdt, D, C, t, L)
            dxdt[1] = -L[1]*C[1]
        end
        
    @inline function A!(a, D, C, t, L)
        a[1] = L[1]*D[1]
    end
    
    model = MassConversion.MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mat, θ, A!, F!)
    
    println("Number of threads ", Threads.nthreads())
    par_run_sim(model, 10_000, 1_000)
end;

reps = Meta.parse(string(ARGS[1]))
blocksize = Meta.parse(string(ARGS[2]))

main(ARGS[1], ARGS[2])