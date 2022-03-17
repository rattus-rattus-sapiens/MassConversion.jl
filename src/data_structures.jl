struct SSAOutput{T}
    trajectories_raw::Array{T, 3}
    mean_total::Array{T, 2}
    var_total::Array{T, 2}
    t_span::Vector{Float64}
    function SSAOutput{T}(ns::Int64, nt::Int64, nr::Int64) where T <: Real
        new(
            Array{T, 3}(undef, (ns, nt, nr)), 
            Array{T, 2}(undef, (ns, nt)),
            Array{T, 2}(undef, (ns, nt)), 
            zeros(Float64, nr)
        )
    end
end

struct MCMOutput{T}
    trajectories_discrete::Array{T, 3}
    trajectories_continuum::Array{T, 3}
    mean_discrete::Array{T, 2}
    mean_continuum::Array{T, 2}
    mean_total::Array{T, 2}
    var_total::Array{T, 2}
    t_span::Vector{Float64}
    function MCMOutput{T}(ns::Int64, nt::Int64, nr::Int64) where T <: Real
        new(
            Array{T, 3}(undef, (ns, nt, nr)),
            Array{T, 3}(undef, (ns, nt, nr)),
            Array{T, 2}(undef, (ns, nt)),
            Array{T, 2}(undef, (ns, nt)),
            Array{T, 2}(undef, (ns, nt)),
            Array{T, 2}(undef, (ns, nt)),
            zeros(Float64, nr)
        )
    end
end

function calc_mean(t::Array{T, 3}) where T <: Real
    ns, nt, nr = size(t)
    mn = Array{T, 2}(undef, (ns, nt))
    for j=1:nt, i=1:ns
        mn[i, j] = sum(t[i,j,:])
    end
    mn ./= nr
end

function calc_mean_total(mn1::Array{T, 2}, mn2::Array{T, 2}) where T <: Real
    ns, _ = size(mn1)
    mnt = similar(mn1)
    for i=1:ns
        @. mnt[i, :] = mn1[i, :] + mn2[i, :]
    end
    mn ./= 2
end

function calc_var(t::Array{T, 3}, mn::Array{T, 2}) where T <: Real
    ns, nt, nr = size(t)
    var = Array{T, 2}(0, (ns, nt))
    for j=1:ns, i=1:nt
        @. var[i, j] += (mn[i, j] - t[i, j, :])^2
    end
    var ./= nr - 1
end