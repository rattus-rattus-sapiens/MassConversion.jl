abstract type Output end

struct SSAOutput <: Output
    trajectories_discrete::Array
    tf::Float64
    dt::Float64
    IC::Vector{Int64}
    λ::Vector{Float64} 
    R::Matrix
    rn::Int64
    description::String
    function SSAOutput(t::Array{T, 3}, tf, dt, IC, λ, R, rn, desc::String="") where T <: Real
        t = permutedims(t, [2, 1, 3])
        new(t, tf, dt, IC, λ, R, rn, desc)
    end
end

struct MCMOutput <: Output
    μd::Array{Float64}
    μc::Array{Float64}
    σd::Array{Float64}
    σc::Array{Float64}
    μt::Array{Float64}
    tspan::Vector{Float64}
    tf::Float64
    dt::Float64
    IC::Vector{Float64}
    λ::Vector{Float64}
    R::Matrix{Int64}
    θ::Vector{Tuple{Float64, Float64}}
    rn::Int64
    description::String
    function MCMOutput(md, mc, Sd, Sc, tf, dt, IC, λ, R, θ, rn, desc::String="")
        md = permutedims(md, [2 1])
        mc = permutedims(mc, [2 1])
        mt = md .+ mc
        Sd = permutedims(Sd, [2 1]) ./ (rn - 1)
        Sc = permutedims(Sc, [2 1]) ./ (rn - 1)
        tspan = Vector(0:dt:tf)
        new(md, mc, Sd, Sc, mt, tspan, tf, dt, IC, λ, R, θ, rn, desc)
    end
end

"""
    to_dict(::O) where O <: Output
"""
function to_dict(dat::O) where {O<:Output}
    return Dict(key => getfield(dat, key) for key ∈ fieldnames(O))
end