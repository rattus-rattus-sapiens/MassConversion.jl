using MassConversion
using Plots

casename = "alt-exp";
filename = "fixed-propensities";
data = load_data(casename, filename);

p = plot_init()
for (id, pair) in enumerate(data)
    plot_rel_err(p, pair[1], pair[2], stride=1, label="")
end
p
#plot_save("rel-errs-1e5x16")

p = plot_init()
plot_total(p, data[1][1], "Total")
plot_both(p, data[1][1])
plot_total(p, data[1][2], "SSA")
hline!([300], color="black", lw=0.5, label="")
#plot_save("profile-1e5x16")

#plot_save()