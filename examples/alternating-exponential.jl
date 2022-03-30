using Revise
using MassConversion
using Plots
using BenchmarkTools

function main()
    tf = 10;
    dt = 1e-2;
    IC = [0, 100];
    λ = [1e0, 1.0];
    R = reshape([-1 0], 2, 1);
    θ = [(50, 50)];

    function F!(dxdt, S, t, L)
        dxdt[2] = -L[1]*S[2]
    end;

    function A!(a, S, t, L)
        a[1] = L[1]*S[1]
    end;

    rn = 100;

    return run_mcm(tf, dt, IC, λ, R, θ, F!, A!, rn)
end

@time rec = main()
md = MassConversion.calc_mean(rec.trajectories_discrete)
mc = MassConversion.calc_mean(rec.trajectories_continuum)

plot(md')
plot!(mc')
plot!(md' .+ mc')