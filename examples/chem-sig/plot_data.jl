using Revise
using Plots
using MassConversion
gr()

dir = "alter";

ssa = load_raw(joinpath("/home/jkynaston/git/MassConversion.jl/examples/alt-log/dat/", dir*"-ssa"));
mcm = load_raw(joinpath("/home/jkynaston/git/MassConversion.jl/examples/alt-log/dat/", dir*"-mcm"));

ssa_full = sum(ssa);
mcm_full = sum(mcm);

p = plot(dpi=600)
plot(p,mcm_full)

err_full = RelativeError(mcm_full, ssa_full)
err = RelativeError(mcm, ssa)

plot(err)
plot(err_full)

t = ode_full.t_range;
ode_dat = ode_full.C;
mcm_dat = mcm_full.D .+ mcm_full.C;
tru_dat = f.(t);

rel_err_mcm = zeros(t.len, length(mcm))

rel_err_ode = ((ode_dat .- tru_dat) ./ tru_dat);
rel_err_mcm = ((mcm_dat .- tru_dat) ./ tru_dat);

for i in 1:length(mcm)
    println(i)
    rel_err_mcm[:,i] = (((mcm[i].C .+ mcm[i].D) .- tru_dat) ./ tru_dat);
end

p = plot(t, rel_err_ode;label="ode",dpi=600);
plot!(p, t[1:100:end], rel_err_mcm[1:100:end,:];label="",dpi=600)
vline!([1.38])