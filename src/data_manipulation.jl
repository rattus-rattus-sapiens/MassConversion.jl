function calc_mean(t::Array{T, 3}) where T <: Real
    nt, ns, nr = size(t)
    mn = zeros(T, nt, ns)
    @inbounds for k=1:nr, j=1:ns, i=1:nt
        mn[i, j] += t[i,j,k]
    end
    mn ./= nr
end

function calc_var(t::Array{T, 3}, mn::Array{T, 2}) where T <: Real
    nt, ns, nr = size(t)
    var = zeros(T, nt, ns)
    @inbounds for k=1:nr, j=1:ns, i=1:nt
        var[i, j] += (mn[i, j] - t[i, j, k])^2
    end
    var ./= nr - 1
end

function to_dict(dat::O) where O <: Output
    return Dict(key=>getfield(dat, key) for key ∈ fieldnames(O))
end