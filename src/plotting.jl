function plot_total(p::P, data::O, label::String="Total") where {P<:Plots.Plot,O<:Output}
    tspan = get_tspan(data)
    mt, = calc_mean(data)
    plot!(p, tspan, mt, label=label)
end

function plot_both(p::P, data::MCMOutput) where {P<:Plots.Plot}
    mt, md, mc = calc_mean(data)
    tspan = get_tspan(data)
    plot!(p, tspan, md, label="Discrete")
    plot!(p, tspan, mc, label="Continuum")
end

function plot_init(size::Tuple{Int64,Int64}=(600, 400))
    pgfplotsx()
    default(legendfontsize=12, tickfontsize=10)
    p = plot(dpi=300, size=size)
    return p
end

function plot_rel_err(rec1::MCMOutput, rec2::MCMOutput)
    mt1 ,= calc_mean(rec1)
    mt2 ,= calc_mean(rec2)
    err = zeros(size(mt1))
    @. err = (mt1 - mt2) / mt2
    tspan = 0:rec1.dt:rec1.tf
    plot!(tspan, err, label="Relative Err")
end

function plot_save(filename::String="")
    if filename == ""
        filename = "examples/plot/" * string(floor(now(), Dates.Second))
    end
    savefig(filename)
    println("File saved at " * filename)
end