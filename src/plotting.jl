function plot_total(data::MCMOutput)
    md = calc_mean(data.trajectories_discrete)
    mc = calc_mean(data.trajectories_continuum)
    tspan = 0:data.dt:data.tf
    plot!(tspan, md .+ mc, label="Total")
end

function plot_total(data::SSAOutput)
    md = calc_mean(data.trajectories_discrete)
    tspan = 0:data.dt:data.tf
    plot!(tspan,md)
end

function plot_both(data::MCMOutput)
    md = calc_mean(data.trajectories_discrete)
    mc = calc_mean(data.trajectories_continuum)
    tspan = 0:data.dt:data.tf
    plot!(tspan, md, label="Discrete")
    plot!(tspan, mc, label="Continuum")
end

function plot_init(size::Tuple{Int64, Int64}=(600,400))
    default(legendfontsize = 12, tickfontsize = 10)
    p = plot(dpi=300, size=size)
    return p
end