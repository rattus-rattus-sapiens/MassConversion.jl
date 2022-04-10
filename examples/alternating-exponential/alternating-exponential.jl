using MassConversion
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--repno", "-r"
            help = "Number of repeats to simulate"
            arg_type = Int
        "--blocksize", "-b"
            help = "Number of repeats per save datafile"
            arg_type = Int
    end
    return parse_args(s)
end

function main(repno, blocksize)
    t_span = 0.0:5e-4:10
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
    par_run_sim(model, repno, blocksize)
end;

parsed_args = parse_commandline()
repno = parsed_args["repno"]
blocksize = parsed_args["blocksize"]

@elapsed main(1_000_000, 10_000)