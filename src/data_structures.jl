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
    params::Dict
    description::String
    function MCMOutput(t::Array{T, 3}, params::Dict) where T <: Real
        ns::Int64 = size(t, 1)/2;
        t1 = t[1:ns, :, :]
        t2 = t[ns+1:end, :, :]
        new(t1, t2, params, "")
    end
end