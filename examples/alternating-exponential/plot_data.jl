using Revise
using Plots
using MassConversion

a, b = load_raw("/home/jkynaston/git/MassConversion.jl/examples/alternating-exponential/dat/Tue-12-00:05:36")

data_mcm_full = sum(a)
data_ssa_full = sum(b)

plot(data_mcm_full.D)
plot!(data_mcm_full.C)
plot!(data_mcm_full.D .+ data_mcm_full.C)

plot!(data_ssa_full.D)

gt = data_ssa_full.D
pt = data_mcm_full.D .+ data_mcm_full.C 

rel_err = (pt .- gt) ./ gt

plot(rel_err)