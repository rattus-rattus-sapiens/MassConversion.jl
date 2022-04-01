function calc_mean(rec::MCMOutput)
    md = dropdims(mean(rec.trajectories_discrete, dims=3), dims=3)
    mc = dropdims(mean(rec.trajectories_continuum, dims=3), dims=3)
    return md .+ mc, md, mc
end

function calc_mean(rec::SSAOutput)
    return dropdims(mean(rec.trajectories_discrete, dims=3))
end

function get_tspan(rec::O) where {O<:Output}
    return 0:rec.dt:rec.tf
end

function calc_var(t::Array{T,3}, mn::Array{T,2}) where {T<:Real}
    nt, ns, nr = size(t)
    var = zeros(T, nt, ns)
    @inbounds for k = 1:nr, j = 1:ns, i = 1:nt
        var[i, j] += (mn[i, j] - t[i, j, k])^2
    end
    var ./= nr - 1
end

function to_dict(dat::O) where {O<:Output}
    return Dict(key => getfield(dat, key) for key ∈ fieldnames(O))
end