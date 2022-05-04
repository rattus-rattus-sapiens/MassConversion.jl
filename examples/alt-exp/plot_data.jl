using Revise
using Plots
using MassConversion

ode_dir = "yes-thresh-ode";
mcm_dir = "yes-thresh-mcm";

ode = load_raw(joinpath("/home/jkynaston/git/MassConversion.jl/examples/alt-exp/dat/", ode_dir));
mcm = load_raw(joinpath("/home/jkynaston/git/MassConversion.jl/examples/alt-exp/dat/", mcm_dir));

ode_full = sum(ode);
mcm_full = sum(mcm);

plot(ode_full)
plot(mcm_full)

function f(t)
    if t ≤ 0
        1000.0
    elseif 0 < t ≤ 5.5
        1000.0*exp(-t)
    elseif t ≥ 15/2
        1082.0849986238989 * exp(-t+7.5)
    else
        200*t - 417.9150013761012
    end
end;

t = ode_full.t_range;
ode_dat = ode_full.C;
mcm_dat = mcm_full.D .+ mcm_full.C;
tru_dat = f.(t);

rel_err_ode = (ode_dat .- tru_dat) ./ tru_dat;
rel_err_mcm = (mcm_dat .- tru_dat) ./ tru_dat;

n1 = sum(abs.(rel_err_mcm))
n2 = sum(abs.(rel_err_ode))

plot(rel_err_ode)
plot!(rel_err_mcm)