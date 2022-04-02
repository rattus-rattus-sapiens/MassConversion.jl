using Revise
using MassConversion
using Plots
using Base.Threads
using JLD2
using BenchmarkTools

Threads.nthreads()

function main(is_ssa::Bool)
    tf = 10.0
    dt = 5e-4
    λ = [1e0, 2e2, 1e1]
    R = [
        -1 1
        0 0
    ]
    if is_ssa
        IC = [1000, 0]
        θ = [(Inf, Inf)]
    else
        IC = [0, 1000]
        θ = [(300, 300)]
    end

    ts = (2.5, 7.5)

    function F!(dxdt, S, t, L)
        if ts[1] < t < ts[2]
            dxdt[2] = 0
        else
            dxdt[2] = -L[1] * S[2]
        end
    end

    function A!(a, S, t, L)
        if ts[1] < t < ts[2]
            a[1] = 0.0
            a[2] = L[2]
        else
            a[1] = L[1] * S[1]
            a[2] = 0.0
        end
    end

    rn = 1000000

    return run_mcm(tf, dt, IC, λ, R, θ, F!, A!, rn)
end;

O = MassConversion.MCMOutput

data = Vector{Tuple{O,O}}()

Threads.@threads for i = 1:16
    @time push!(data, (main(false), main(true)))
    println("Repeat " * string(i) * " Complete")
end

save_mcm(data)

p = plot_init()
for (id, pair) in enumerate(data)
    plot_rel_err(p, pair[1], pair[2], stride=50, label="")
end
p
plot_save("altexp/rel-errs-SSAvODE")

p = plot_init()
plot_total(p, data[1][1], "Total")
plot_both(p, data[1][1])
hline!([250], color="black", lw=0.5, label="")
plot_save("altexp/profile-ODE")

#plot_save()