abstract type Output end

struct SSAOutput <: Output
    trajectories_discrete::Array
    function SSAOutput(t::Array{T, 3}) where T <: Real
        new(t)
    end
end

struct MCMOutput <: Output
    trajectories_discrete::Array
    trajectories_continuum::Array
    tf::Float64
    dt::Float64
    IC::Vector{Float64}
    λ::Vector{Float64}
    R::Matrix{Int64}
    θ::Vector{Tuple}
    rn::Int64
    description::String
    function MCMOutput(t::Array{T, 3}, tf, dt, IC, λ, R, θ, rn) where T <: Real
        ns = floor(Int64, size(t, 1)/2);
        t1 = t[1:ns, :, :]
        t2 = t[ns+1:end, :, :]
        new(t1, t2, tf, dt, IC, λ, R, θ, rn, "")
    end
end