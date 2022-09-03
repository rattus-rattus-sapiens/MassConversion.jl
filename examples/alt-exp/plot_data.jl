using Revise
using Plots
using MassConversion
pgfplotsx()

dir = "alter";

ode = load_raw(joinpath("/home/jkynaston/git/MassConversion.jl/examples/alt-exp/dat/", dir*"-ode"));
mcm = load_raw(joinpath("/home/jkynaston/git/MassConversion.jl/examples/alt-exp/dat/", dir*"-mcm"));

ode_full = sum(ode);
mcm_full = sum(mcm);

p = plot(dpi=300)
plot(p,mcm_full)
plot!(p, size=(335,200), legend=false, tickfontsize=10)
plot!(legend=false, tickfontsize=10)

function f(t)
    if t ≤ 0
        1000.0
    elseif 0 < t < 4
        1000.0*exp(-t)
    else
        200*(t-4) + 18.315638888734178
    end
end;

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