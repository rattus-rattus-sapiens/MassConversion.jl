using Revise
using Plots
using MassConversion

a, b = load_raw("/home/jkynaston/git/MassConversion.jl/examples/alternating-exponential/dat/Tue-12-19:08:09")

data_mcm_full = sum(a)
data_ssa_full = sum(b)

p = plot(data_mcm_full)
plot(p, data_ssa_full)

gt = data_ssa_full.D
pt = data_mcm_full.D .+ data_mcm_full.C 

rel_err = (pt .- gt) ./ gt

plot(rel_err)