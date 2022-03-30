function calc_mean(t::Array{T, 3}) where T <: Real
    ns, nt, nr = size(t)
    mn = Array{T, 2}(undef, (ns, nt))
    for j=1:nt, i=1:ns
        mn[i, j] = sum(t[i,j,:])
    end
    mn ./= nr
end

function calc_var(t::Array{T, 3}, mn::Array{T, 2}) where T <: Real
    ns, nt, nr = size(t)
    var = zeros(T, ns, nt)
    for i=1:ns, j=1:nt, k=1:nr 
        var[i, j] += (mn[i, j] - t[i, j, k])^2
    end
    var ./= nr - 1
end