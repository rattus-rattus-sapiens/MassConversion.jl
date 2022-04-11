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
    D_mcm = [0]
    C_mcm = [1000]
    λ_reac = [1e0, 0.0]
    λ_tran = 1.0
    R_mcm = [
        -1 0
        0 0
    ]
    θ = [(250, 250)]

    @inline function dxdt_mcm(dxdt, D, C, t, L)
        dxdt[1] = -L[1] * C[1]
    end

    @inline function prop_mcm(a, D, C, t, L)
        a[1] = L[1] * D[1]
    end

    mcm_model = MCMmodel(t_span, D_mcm, C_mcm, λ_reac, λ_tran, R_mcm, θ, prop_mcm, dxdt_mcm)

    D_ssa = [1000]
    R_ssa = Matrix{Int64}(undef, 1, 2)
    R_ssa[1] = -1
    R_ssa[2] = 0

    @inline function prop_ssa(a, D, t, L)
        a[1] = L[1] * D[1]
    end

    ssa_model = SSAmodel(t_span, D_ssa, λ_reac, R_ssa, prop_ssa)

    println("Number of threads ", Threads.nthreads())
    par_run_sim(mcm_model, ssa_model, repno, blocksize)
end;

parsed_args = parse_commandline()
repno = parsed_args["repno"]
blocksize = parsed_args["blocksize"]

main(repno, blocksize)