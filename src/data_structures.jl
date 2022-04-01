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
    trajectories_discrete::Array
    trajectories_continuum::Array
    md::Array
    mc::Array
    tspan::Vector{Float64}
    tf::Float64
    dt::Float64
    IC::Vector{Float64}
    λ::Vector{Float64}
    R::Matrix{Int64}
    θ::Vector{Tuple}
    rn::Int64
    description::String
    function MCMOutput(t::Array{T, 3}, tf, dt, IC, λ, R, θ, rn, desc::String="") where T <: Real
        ns = floor(Int64, size(t, 1)/2);
        t1 = permutedims(t[1:ns, :, :], [2, 1, 3])
        t2 = permutedims(t[ns+1:end, :, :], [2, 1, 3])
        md = dropdims(mean(t1, dims=3), dims=3)
        mc = dropdims(mean(t2, dims=3), dims=3)
        tspan = Vector(0:dt:tf)
        new(t1, t2, md, mc, tspan, tf, dt, IC, λ, R, θ, rn, desc)
    end
end