function plot_total(p::P, data::O, label::String="Total") where {P<:Plots.Plot,O<:Output}
    mt = data.μd .+ data.μc
    plot!(p, data.tspan, mt, label=label)
end

function plot_both(p::P, data::MCMOutput) where {P<:Plots.Plot}
    plot!(p, data.tspan, data.μd , label="Discrete")
    plot!(p, data.tspan, data.μc , label="Continuum")
end

function plot_init(size::Tuple{Int64,Int64}=(600, 400))
    pgfplotsx()
    default(legendfontsize=12, tickfontsize=10)
    p = plot(dpi=300, size=size)
    return p
end

function plot_rel_err(p::Plots.Plot, rec1::MCMOutput, rec2::MCMOutput; stride::Int64=1, label::String="Relative Err")
    mt1 = rec1.μd[1:stride:end] .+ rec1.μc[1:stride:end]
    mt2 = rec2.μd[1:stride:end] .+ rec2.μc[1:stride:end]
    err = similar(mt1)
    @. err = (mt1 - mt2) / mt2
    plot!(p, rec1.tspan[1:stride:end], err, label=label)
end

function plot_abs_err(rec1::MCMOutput, rec2::MCMOutput)
    mt1 = rec1.μd .+ rec1.μc
    mt2 = rec2.μd .+ rec2.μc
    plot!(rec1.tspan, mt1 .- mt2, label="Absolute error")
end

function plot_save(filename::String="")
    if filename == ""
        filename = "examples/plot/" * string(floor(now(), Dates.Second))
    else
        filename = "examples/plot/" * filename
    end
    savefig(filename)
    println("File saved at " * filename)
end