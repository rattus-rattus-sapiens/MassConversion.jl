using MassConversion
using ArgParse
using Dates

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--repno", "-r"
        help = "Number of repeats to simulate"
        arg_type = Int
        "--blocksize", "-b"
        help = "Number of repeats per save datafile"
        arg_type = Int
        "--dir_name", "-d"
        help = "Name of saved data directory"
        arg_type = String
        default = Dates.format(now(), "e-dd-HH:MM:SS")
    end
    return parse_args(s)
end

function main(repno, blocksize, dir_name)
    t_span = 0.0:1e-2:10
    D_mcm = [0]
    C_mcm = [1000]
    λ_reac = [1e0, 2e2]
    λ_tran = 1.0
    R_mcm = [
        -1 1
        0 0
    ]
    θ = [(0, 0)]

    @inline function dxdt_mcm(dxdt, D, C, t, L)
        if 2.5 ≤ t ≤ 7.5
            dxdt[1] = 0
        else
            dxdt[1] = -L[1] * C[1]
        end
    end

    @inline function prop_mcm(a, D, C, t, L)
        if 2.5 ≤ t ≤ 7.5
            a[1] = 0.0
            a[2] = L[2]
        else
            a[1] = L[1] * D[1]
            a[2] = 0.0
        end
    end

    mcm_model = MCMmodel(t_span, D_mcm, C_mcm, λ_reac, λ_tran, R_mcm, θ, prop_mcm, dxdt_mcm)

    D_ssa = [1000]
    R_ssa = Matrix{Int64}(undef, 1, 2)
    R_ssa[1] = -1
    R_ssa[2] = 1

    @inline function prop_ssa(a, D, t, L)
        if 2.5 ≤ t ≤ 7.5
            a[1] = 0.0
            a[2] = L[2]
        else
            a[1] = L[1] * D[1]
            a[2] = 0
        end
    end

    ssa_model = SSAmodel(t_span, D_ssa, λ_reac, R_ssa, prop_ssa)

    println("Number of threads ", Threads.nthreads())
    par_run_sim(mcm_model, ssa_model, repno, blocksize; dir_name=dir_name)
end

args = parse_commandline()
main(args["repno"], args["blocksize"], args["dir_name"])