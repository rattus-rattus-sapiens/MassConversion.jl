println("Initialising...")

using MassConversion
using Base.Threads
using DifferentialEquations

println("Number of threads ", Threads.nthreads())

function main(is_ssa::Bool, rn::Int64)
    tf = 10.0
    λ = [1e0, 2e2, 1e2]
    dt = 5e-4
    R = [
        -1 1
        0 0
    ]
    if is_ssa
        IC = [1000, 0]
        θ = [(Inf, Inf)]
    else
        IC = [0, 1000]
        θ = [(0, 0)]
    end

    ts = (2.5, 7.5)

    function F!(dxdt, S, t, L)
        if is_ssa
            dxdt[2] = 0.0
        else
            if ts[1] ≤ t ≤ ts[2]
                dxdt[2] = L[2] 
            else
                dxdt[2] = -L[1] * S[2]
            end
        end
    end

    function A!(a, S, t, L)
        if is_ssa
            if ts[1] ≤ t ≤ ts[2]
                a[1] = 0.0
                a[2] = L[2]
            else
                a[1] = L[1] * S[1]
                a[2] = 0.0
            end
        else
            a[1] = 0.0
            a[2] = 0.0
        end
        return nothing
    end

    return run_mcm(tf, dt, IC, λ, R, θ, F!, A!, rn)
end;

O = MassConversion.MCMOutput;
data = Vector{Tuple{O,O}}();

Threads.@threads for i = 1:16
    time = @elapsed push!(data, (main(false, 100000), main(true, 100000)))
    println("Repeat " * string(i) * " complete. Elapsed: " * string(time) * " seconds.")
end

#! IO CONFIG
casename = "alt-exp";
datname = "fixed-propensities";
save_data(data, casename, datname);