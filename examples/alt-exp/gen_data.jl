using MassConversion

function gen_ode_data(repno, blocksize, dir_name)
    println("Beginning ODE simulation")
    t_span = 0.0:1e-4:8.0
    D_init = [0]
    C_init = [1000]
    λ_reac = [1e0, 2e2]
    λ_tran = 0.0
    R_mcm = zeros(Int64,2,2)
    theta = [(0,0)]
    @inline function dxdt_mcm(k, t, C)
        @inbounds if t ≥ 4.0
            k[1] = λ_reac[2]
        else
            k[1] = C[1] * -λ_reac[1]
        end
        return nothing
    end
    @inline function prop_mcm(a, D, C, t, L)
        return nothing
    end
    model = MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mcm, theta, prop_mcm, dxdt_mcm)
    par_run_sim(model, 1, 1; dir_name=dir_name)
end

function gen_mcm_data(repno, blocksize, dir_name)
    println("Beginning ODE simulation")
    t_span = 0.0:1e-4:8.0
    D_init = [0]
    C_init = [1000]
    λ_reac = [1e0, 2e2]
    λ_tran = 1.0e1
    R_mcm = [
        -1 1
        0 0
        ]
    theta = [(650,700)]
    @inline function f(k, t, C)
        @inbounds if t ≥ 4.0
            k[1] = 0.0
        else
            k[1] = C[1] * -λ_reac[1]
        end
        return nothing
    end
    @inline function prop_mcm(a, D, C, t, L)
        @inbounds if t ≥ 4.0
            a[1] = 0
            a[2] = L[2]
        else
            a[1] = L[1]*D[1]
            a[2] = 0
        end
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