using Revise
using Plots
using MassConversion
pgfplotsx()

dir = "meeting";

ode = load_raw(joinpath("/home/jkynaston/git/MassConversion.jl/examples/alt-exp/dat/", dir*"-ode"));
mcm = load_raw(joinpath("/home/jkynaston/git/MassConversion.jl/examples/alt-exp/dat/", dir*"-mcm"));

ode_full = sum(ode);
mcm_full = sum(mcm);

p = plot(dpi=600)
plot(p,mcm_full)

function f(t)
    if t ≤ 0
        1000.0
    elseif 0 < t ≤ 2.5
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

rel_err_ode = ((ode_dat .- tru_dat) ./ tru_dat);
rel_err_mcm = ((mcm_dat .- tru_dat) ./ tru_dat);

n1 = sum(abs.(rel_err_mcm))
n2 = sum(abs.(rel_err_ode))

p = plot(t, rel_err_ode;label="ode",dpi=600);
plot!(p, t, rel_err_mcm;label="mcm");
p
vline!([1.38])