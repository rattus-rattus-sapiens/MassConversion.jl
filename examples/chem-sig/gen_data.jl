using MassConversion

function gen_mcm_data(repno, blocksize, dir_name)
    println("Beginning ODE simulation")
    t_span = 0.0:1e-1:1000.0
    D_init = [10, 100, 0]
    C_init = [0, 0, 0]
    λ_reac = [1e2, 1e3, 1e-2, 9900.0, 9900.0, 9900.0, 5e2, 1e2]
    λ_tran = 1.0e0
    R_mcm = [
        0 0 0 0 0 0 1 -1
        0 1 -1 0 0 0 0 0
        1 -1 0 -1 0 -1 0 0
        0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0
        0 0 0 0 -1 0 0 0
    ]
    theta = [(Inf, Inf), (200, 250), (Inf, Inf)]
    @inline function f(k, t, C)
        k[1] = -λ_reac[8] * C[1]
        k[2] = λ_reac[2] * C[3] - λ_reac[3] * C[2]
        k[3] = -λ_reac[2] * C[3] - λ_reac[4] * C[1] * C[3]
    end
    @inline function prop_mcm(a, D, C, t, L)
        a[1] = λ_reac[1]
        a[2] = λ_reac[2] * D[3]
        a[3] = λ_reac[3] * D[2]
        a[4] = λ_reac[4] * C[1] * D[3]
        a[5] = λ_reac[5] * D[1] * C[3]
        a[6] = λ_reac[6] * D[1] * D[3]
        a[7] = λ_reac[7]
        a[8] = λ_reac[8] * D[1]
    end
    model = MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mcm, theta, prop_mcm, f)
    par_run_sim(model, repno, blocksize; dir_name=dir_name)
end

function gen_ssa_data(repno, blocksize, dir_name)
    t_span = 0.0:1e-1:1000.0
    D_init = [10, 100, 0]
    λ_reac = [1e2, 1e3, 1e-2, 9900.0, 5e2, 1e2]
    R_ssa = [
        0 0 0 0 1 -1
        0 1 -1 0 0 0
        1 -1 0 -1 0 0
    ]
    @inline function prop_ssa(a, D, t, L)
        a[1] = λ_reac[1]
        a[2] = λ_reac[2] * D[3]
        a[3] = λ_reac[3] * D[2]
        a[4] = λ_reac[4] * D[1] * D[3]
        a[5] = λ_reac[5]
        a[6] = λ_reac[6] * D[1]
    end
    model = SSAmodel(t_span, D_init, λ_reac, R_ssa, prop_ssa)
    par_run_sim(model, repno, blocksize; dir_name=dir_name)
end

function gen_ode_data(dir_name)
    println("Beginning ODE simulation")
    t_span = 0.0:1e-1:1000.0
    D_init = [0, 0, 0]
    C_init = [10, 100, 0]
    λ_reac = [1e2, 1e3, 1e-2, 9900, 5e2, 1e2, 0e0, 0e0]
    λ_tran = 0.0e0
    R_mcm = [
        0 0 0 0 0 0 1 -1
        0 1 -1 0 0 0 0 0
        1 -1 0 -1 0 -1 0 0
        0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0
        0 0 0 0 -1 0 0 0
    ]
    theta = [(Inf, Inf), (200, 250), (Inf, Inf)]
    @inline function f(k, t, C)
        a, b, c = C
        k[1] = λ_reac[5] - λ_reac[6] * a
        k[2] = λ_reac[2] * c - λ_reac[3] * b
        k[3] = λ_reac[1] - λ_reac[2] * c - λ_reac[4] * a * c
    end
    @inline function prop_mcm(a, D, C, t, L)
        a[1] = 0#λ_reac[1]
        a[2] = 0#λ_reac[2]*D[3]
        a[3] = 0#λ_reac[3]*D[2]
        a[4] = 0#λ_reac[4]*C[1]*D[3]
        a[5] = 0#λ_reac[5]*D[1]*C[3]
        a[6] = 0#λ_reac[6]*D[1]*D[3]
        a[7] = 0#λ_reac[7]
        a[8] = 0#λ_reac[8]*D[1]
    end
    model = MCMmodel(t_span, D_init, C_init, λ_reac, λ_tran, R_mcm, theta, prop_mcm, f)
    par_run_sim(model, 1, 1; dir_name=dir_name)
end


args = parse_cmd()
repno = args["repno"]
bs = args["blocksize"]
dir = args["dir_name"]
gen_ssa_data(repno, bs, dir * "-ssa")
gen_mcm_data(repno, bs, dir * "-mcm")
gen_ode_data(dir * "-ode")