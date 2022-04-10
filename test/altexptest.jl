println("Initialising...")

using MassConversion
using Base.Threads

println("Number of threads ", Threads.nthreads())

function main(is_ssa::Bool, reac_num::Int64, pid::Int64)
    t_fin = 10.0
    λ_reac = [1e0, 2e2, 1e2]
    t_step = 5e-4
    reac_mat = [
        -1 1
        0 0
    ]
    if is_ssa
        C_init = [1000, 0]
        θ_list = [(Inf, Inf)]
    else
        C_init = [0, 1000]
        θ_list = [(0, 0)]
    end

    ts = (2.5, 7.5)

    function F!(dxdt, S, t, λ_reac)
        if is_ssa
            dxdt[2] = 0.0
        else
            if ts[1] ≤ t ≤ ts[2]
                dxdt[2] = λ_reac[2] 
            else
                dxdt[2] = -λ_reac[1] * S[2]
            end
        end
    end

    function A!(a, S, t, λ_reac)
        if is_ssa
            if ts[1] ≤ t ≤ ts[2]
                a[1] = 0.0
                a[2] = λ_reac[2]
            else
                a[1] = λ_reac[1] * S[1]
                a[2] = 0.0
            end
        else
            a[1] = 0.0
            a[2] = 0.0
        end
        return nothing
    end

    return run_mcm(t_fin, t_step, C_init, λ_reac, reac_mat, θ_list, F!, A!, reac_num)
end;

O = MassConversion.MCMOutput;
data = Vector{Tuple{O,O}}();

Threads.@threads for i = 1:16
    time = @elapsed push!(data, (main(false, 100), main(true, 100)))
    println("Repeat " * string(i) * " complete. Elapsed: " * string(time) * " seconds.")
end

#! IO CONFIG
casename = "alt-exp";
datname = "fixed-propensities";
save_data(data, casename, datname);