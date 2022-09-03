using MassConversion

function gen_mcm_data(repno, blocksize, dir_name)
    println("Beginning ODE simulation")
    t_span = 0.0:1e-2:80.0
    D_init = [60]
    C_init = [0]
    λ_reac = [1e-3, 1e-3, 1e-3, 0.6]
    λ_tran = 1.0e0
    R_mcm = [
        -1 -1 -1 1
        0 0 0 0
        ]
    theta = [(300,300)]
    @inline function f(k, t, C)
        @inbounds if t ≤ 20 || 35 < t ≤ 55 || 70 < t
            k[1] = λ_reac[4]*C[1] - λ_reac[1]*(C[1]^2)
        else
            k[1] = -λ_reac[1]*(C[1]^2)
        end
        return nothing
    end
    @inline function prop_mcm(a, D, C, t, L)
        @inbounds if t ≤ 20 || 35 < t ≤ 55 || 70 < t
            a[1] = L[1] * D[1] * (D[1] - 1)
            a[2] = L[2] * D[1] * C[1]
            a[3] = L[3] * D[1] * C[1]
            a[4] = L[4] * D[1]
        else
            a[1] = L[1] * D[1] * (D[1] - 1)
            a[2] = L[2] * D[1] * C[1]
            a[3] = L[3] * D[1] * C[1]
            a[4] = 0.0
        end
    end
    model = MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mcm, theta, prop_mcm, f)
    par_run_sim(model, repno, blocksize; dir_name=dir_name)
end

function gen_ode_data(dir_name)
    println("Beginning ODE simulation") 
    t_span = 0.0:1e-2:80.0
    D_init = [0]
    C_init = [60]
    λ_reac = [1e-3, 1e-3, 1e-3, 0.6]
    λ_tran = 1.0e0
    R_mcm = [
        -1 -1 -1 1
        0 0 0 0
        ]
    theta = [(0,0)]
    @inline function f(k, t, C)
        @inbounds if t ≤ 20 || 35 < t ≤ 55 || 70 < t
            k[1] = λ_reac[4]*C[1] - λ_reac[1]*(C[1]^2)
        else
            k[1] = -λ_reac[1]*(C[1]^2)
        end
        return nothing
    end
    @inline function prop_mcm(a, D, C, t, L)

    end
    model = MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mcm, theta, prop_mcm, f)
    par_run_sim(model, 1, 1; dir_name=dir_name)
end

function gen_ssa_data(repno, blocksize, dir_name)
    t_span = 0.0:1e-2:80.0
    D_init = [60]
    λ_reac = [1e-3, 0.6]
    R_ssa = [-1;;1]
    @inline function prop_ssa(a, D, t, L)
        @inbounds if t ≤ 20 || 35 < t ≤ 55 || 70 < t
            a[1] = L[1] * D[1] * (D[1] - 1)
            a[2] = L[2] * D[1]
        else
            a[1] = L[1] * D[1] * (D[1] - 1)
            a[2] = 0.0
        end
        return nothing
    end
    model = SSAmodel(t_span, D_init, λ_reac, R_ssa, prop_ssa)
    par_run_sim(model, repno, blocksize; dir_name=dir_name)
end

args = parse_cmd()
repno = args["repno"]
bs = args["blocksize"]
dir = args["dir_name"]
gen_ssa_data(repno, bs, dir*"-ssa1")
gen_ssa_data(repno, bs, dir*"-ssa2")
gen_mcm_data(repno, bs, dir*"-mcm")
gen_ode_data(dir*"-ode")