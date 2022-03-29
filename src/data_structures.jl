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
        new(t1, t2, mn1', mn2', mn3', vr1')
    end
end

abstract type AbstractUInterval end

struct Interval <: AbstractUInterval
    left::Float64
    right::Float64
    function Interval(int::Vector{T}) where T <: Real
        if length(int) ≠ 2 throw("not a valid interval") end
        if int[1] > int[2] 
            new(NaN, NaN)
        else
            left::Float64 = int[1]
            right::Float64 = int[2]
            new(left,right)
        end
    end
    function Interval(int::Vector{Any})
        if length(int) == 0 
            new(NaN, NaN)
        else
            throw("not a valid interval")
        end
    end
end

function Base.iterate(itvl::Interval) return (itvl, 1) end
function Base.iterate(itvl::Interval, state) return nothing end

function in(num::T, itvl::Interval) where T <: Real
    return itvl.left <= num <= itvl.right
end

function union(i1::Interval, i2::Interval)
    a = i1.left
    b = i1.right
    c = i2.left
    d = i2.right
    if a <= c < d <= b
        return UInterval([a, b])
    elseif c <= a < b <= d
        return UInterval([c, d])
    elseif c <= a < d <= b
        return UInterval([c, b])
    elseif a <= c < b <= d
        return UInterval([a, d])
    elseif a < b == c < d
        return UInterval([a,d])
    elseif c < d == a < b
        return UInterval([c, b])
    elseif isnan(a) || isnan(b)
        return UInterval([i2])
    elseif isnan(c) || isnan(d)
        return UInterval([i1])
    else
        return UInterval([i1, i2])
    end
end

struct UInterval 
    itvls::Vector{Interval}
    function UInterval(itvl::Interval)
        new([itvl])
    end
    function UInterval(itvl::Vector{T}) where T <: Real
        new([Interval(itvl)]) 
    end
    function UInterval(itvls::Vector{Interval})
        new(itvls)
    end
end

function Base.iterate(uitvl::UInterval) iterate(uitvl.itvls) end
function Base.iterate(uitvl::UInterval, state) iterate(uitvl.itvls, state) end

function in(num::T, itvl::UInterval) where T <: Real
    for itvl in itvl.itvls
        if in(num, itvl)
            return true
            break
        end
    end
    return false
end

function union(uivl::UInterval, ivl::Interval)
    if length(uivl.itvls) == 1
        return uivl.itvls[1] ∪ ivl
    elseif length(uivl.itvls) == 2
        a = uivl.itvls[1] ∪ ivl
        b = uivl.itvls[2] ∪ ivl
        c = a ∪ b
        return c
    else
        newivl = UInterval(ivl)
        for i in 1:length(uivl.itvls)
            newivl = newivl ∪ uivl.itvls[i] 
        end
        return newivl
    end
end

function union(ivl::Interval, uivl::UInterval) return union(uivl, ivl) end

function union(uivl1::UInterval, uivl2::UInterval)
    if length(uivl1.itvls) == 1 && length(uivl2.itvls) == 1
        return uivl1.itvls[1] ∪ uivl2.itvls[1]
    else
        newivl = uivl1
        for ivl in uivl2.itvls  
            newivl = newivl ∪ ivl
        end
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
    return mnt
end

function calc_var(t::Array{T, 3}, mn::Array{T, 2}) where T <: Real
    ns, nt, nr = size(t)
    var = zeros(T, ns, nt)
    for i=1:ns, j=1:nt, k=1:nr 
        var[i, j] += (mn[i, j] - t[i, j, k])^2
    end
    var ./= nr - 1
end