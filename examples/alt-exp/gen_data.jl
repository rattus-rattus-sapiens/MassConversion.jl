using MassConversion

function gen_ode_data(repno, blocksize, dir_name)
    println("Beginning ODE simulation")
    t_span = 0.0:1e-3:5.0
    D_init = [0]
    C_init = [1000]
    λ_reac = [1e0, 2e2]
    λ_tran = 0.0
    R_mcm = zeros(Int64,2,2)
    theta = [(0,0)]
    @inline function dxdt_mcm(dxdt, D, C, t, L)
        if 5.5 ≤ t ≤ 7.5
            dxdt[1] = L[2] 
        else
            dxdt[1] = -L[1] * C[1]
        end
    end
    @inline function prop_mcm(a, D, C, t, L)
        a[1] = 0.0
        a[2] = 0.0
    end
    model = MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mcm, theta, prop_mcm, dxdt_mcm)
    par_run_sim(model, 1, 1; dir_name=dir_name)
end

function gen_mcm_data(repno, blocksize, dir_name)
    println("Beginning ODE simulation")
    t_span = 0.0:1e-3:5.0
    D_init = [0]
    C_init = [1000]
    λ_reac = [1e0, 2e2]
    λ_tran = 1.0e+1
    R_mcm = [
        -1 1
        0 0
    ]
    theta = [(200,300)]
    @inline function dxdt_mcm(dxdt, D, C, t, L)
        if 5.5 ≤ t ≤ 7.5
            dxdt[1] = 0.0
        else
            dxdt[1] = -L[1] * C[1]
        end
    end
    @inline function prop_mcm(a, D, C, t, L)
        if 5.5 ≤ t ≤ 7.5
            a[1] = 0.0
            a[2] = L[2]
        else
            a[1] = L[1]*D[1]
            a[2] = 0
        end
    end
    model = MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mcm, theta, prop_mcm, dxdt_mcm)
    par_run_sim(model, repno, blocksize; dir_name=dir_name)
end

args = parse_cmd()
repno = args["repno"]
bs = args["blocksize"]
dir = args["dir_name"]
gen_ode_data(repno, bs, dir*"-ode")
gen_mcm_data(repno, bs, dir*"-mcm")