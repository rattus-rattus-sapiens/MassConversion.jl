struct SSAOutput
    trajectories_discrete::Array
    mean_total::Array
    var_total::Array
    reac_count::Vector{Int64}
    function SSAOutput(t::Array{T, 3}, r::Vector{Int64}) where T <: Real
        mn = calc_mean(t)
        vr = calc_var(t, mn)
        new(t, mn, vr, r)
    end
end

struct MCMOutput
    trajectories_discrete::Array
    trajectories_continuum::Array
    mean_discrete::Array
    mean_continuum::Array
    mean_total::Array
    var_total::Array
    function MCMOutput(t::Array{T, 3}) where T <: Real
        if size(t, 1)/2 ≠ floor(size(t, 1)/2) throw(ArgumentError("inital condition must be 2* number of species")) end
        ns::Int64 = size(t, 1)/2;
        t1 = t[1:ns, :, :]
        t2 = t[ns+1:end, :, :]
        mn1 = calc_mean(t1)
        mn2 = calc_mean(t2)
        mn3 = calc_mean_total(mn1, mn2)
        vr1 = calc_var(t1 .+ t2, mn3)
        new(t1, t2, mn1, mn2, mn3, vr1)
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
    mnt ./= 2
end

function calc_var(t::Array{T, 3}, mn::Array{T, 2}) where T <: Real
    ns, nt, nr = size(t)
    var = zeros(T, ns, nt)
    for i=1:ns, j=1:nt, k=1:nr 
        var[i, j] += (mn[i, j] - t[i, j, k])^2
    end
    var ./= nr - 1
end