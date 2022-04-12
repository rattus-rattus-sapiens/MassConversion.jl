# ! ----- SSA Model Type -----
struct SSAmodel{F,S1,S2} <: AbstractModel
    t_start::Float64
    t_step::Float64
    t_final::Float64
    t_len::Int64
    n_spec::Int64
    D_init::SVector{S1,Int64}
    n_reac::Int64
    R_mat::Matrix{Int64}
    λ_reac::SVector{S2,Float64}
    calc_prop::F

    # Constructor
    function SSAmodel(
        t_span::AbstractRange,
        D_init::AbstractVector{Int64},
        λ_reac::AbstractVector{T1},
        R_mat::Matrix{Int64},
        calc_prop::Function
    ) where {T1<:Real}

        t_start = Float64(t_span.ref)
        t_step = Float64(t_span.step)
        t_len = t_span.len
        t_final = t_start + t_step * (t_len - 1)
        n_spec = length(D_init)
        n_reac = size(R_mat, 2)

        # typecast
        D_init = SVector{n_spec,Int64}(D_init)
        λ_reac = SVector{n_reac,Float64}(λ_reac)

        new{typeof(calc_prop),n_spec,n_reac}(
            t_start, t_step, t_final, t_len, n_spec, D_init, n_reac, R_mat, λ_reac, calc_prop
        )
    end
end

get_t_range(S::SSAmodel) = S.t_start:S.t_step:S.t_final

# ! ----- SSA Raw Data Type -----
mutable struct SSAraw{T}
    D::Array{Int64,2}
    n::Int64
    t_range::T
    # Constructors
    function SSAraw(S::SSAmodel) 
        tr = get_t_range(S)
        new{typeof(tr)}(zeros(Int64, S.n_spec, S.t_len), 0, tr)
    end
end

@inline function record!(raw::SSAraw, D, ti)
    @inbounds @simd for k in 1:length(D)
        raw.D[k, ti] += D[k]
    end
end

# ! ----- SSA Actual Data Type -----
struct SSAdata{T}
    D::Array{Float64,2}
    n::Int64
    t_range::T
    # Constructors
    function SSAdata(raw::SSAraw)
        tr = raw.t_range
        new{typeof(tr)}(raw.D' ./ raw.n, raw.n, raw.t_range)
    end
    SSAdata(D::Array{Float64,2}, n::Int64, tr) = new{typeof(tr)}(D, n, tr)
end

Base.zero(::Type{SSAdata}) = SSAdata(zeros(Float64, 0, 0), 0, 0:0.1:0)
Base.iszero(data::SSAdata) = isempty(data.D)

# ? --- SSA Actual Data functions ---
function +(d1::SSAdata, d2::SSAdata)
    if iszero(d1)
        return d2
    elseif iszero(d2)
        return d1
    else
        n1 = d1.n
        n2 = d2.n
        n3 = n1 + n2
        D = similar(d1.D)
        @. D = (d1.D * n1 + d2.D * n2) / n3
        return SSAdata(D, n3, d1.t_range)
    end
end