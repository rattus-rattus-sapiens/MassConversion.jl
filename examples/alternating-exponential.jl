using Revise
using MassConversion
using Plots
using BenchmarkTools

function main(is_ssa::Bool)
    tf = 10.0
    dt = 1e-2
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
        θ = [(100, 100)]
    end

    ts = (5.0, 7.5)

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

    rn = 1000

    return run_mcm(tf, dt, IC, λ, R, θ, F!, A!, rn)
end;

rec1 = main(false);
rec2 = main(true);

p = plot_init()
plot_total(p, rec1, "Total")
plot_both(p, rec1)
plot_total(p, rec2, "Pure SSA")

p = plot_init()
plot_rel_err(rec1, rec2)
plot_save()