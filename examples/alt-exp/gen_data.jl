using MassConversion

function gen_ode_data(repno, blocksize, dir_name)
    println("Beginning ODE simulation")
    t_span = 0.0:1e-4:5.0
    D_init = [0]
    C_init = [1000]
    λ_reac = [1e0]
    λ_tran = 0.0
    R_mcm = zeros(Int64,2,1)
    theta = [(0,0)]
    @inline function dxdt_mcm(t, C)
        return C .* -λ_reac[1]
    end
    @inline function prop_mcm(a, D, C, L)
        return nothing
    end
    model = MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mcm, theta, prop_mcm, dxdt_mcm)
    par_run_sim(model, 1, 1; dir_name=dir_name)
end

function gen_mcm_data(repno, blocksize, dir_name)
    println("Beginning ODE simulation")
    t_span = 0.0:1e-4:5.0
    D_init = [0]
    C_init = [1000]
    λ_reac = [1e0, 0.0]
    λ_tran = 1.0e1
    R_mcm = [
        -1 1
        0 0
        ]
    theta = [(650,700)]
    @inline function f(t, C)
        return C .* -λ_reac[1]
    end

    @inline function prop_mcm(a, D, C, L)
        a[1] = L[1]*D[1]
        a[2] = 0
    end
    model = MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mcm, theta, prop_mcm, f)
    par_run_sim(model, repno, blocksize; dir_name=dir_name)
end

args = parse_cmd()
repno = args["repno"]
bs = args["blocksize"]
dir = args["dir_name"]
gen_ode_data(repno, bs, dir*"-ode")
gen_mcm_data(repno, bs, dir*"-mcm")