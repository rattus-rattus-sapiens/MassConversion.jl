using Revise
using Plots
using MassConversion

_, ssa = load_raw("/home/jkynaston/git/MassConversion.jl/examples/alt-exp/dat/ssa-solution");
ode, _ = load_raw("/home/jkynaston/git/MassConversion.jl/examples/alt-exp/dat/ode-solution");
mcm, _ = load_raw("/home/jkynaston/git/MassConversion.jl/examples/alt-exp/dat/mcm-solution");

ssa_full = sum(ssa);
ode_full = ode[1];
mcm_full = sum(mcm);

plot(ode_full)
plot(ssa_full)
plot(mcm_full)

err_ode = RelativeError(ode_full, ssa_full)
err_mcm = RelativeError(mcm_full, ssa_full)
p = plot(err_ode)
plot(p, err_mcm)

function f(t)
    if t ≤ 0
        1000.0
    elseif 0 < t ≤ 5.5
        1000.0*exp(-t)
    elseif t > 15/2
        1082.0849986238989 * exp(-t+7.5)
    else
        200*t - 417.9150013761012
    end
end

gt = data_ssa_full.D
pt = data_mcm_full.D .+ data_mcm_full.C 

rel_err = (pt .- gt) ./ gt

plot(rel_err)