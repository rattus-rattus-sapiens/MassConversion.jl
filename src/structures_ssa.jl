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

# ! ----- SSA Raw Data Type -----
mutable struct SSAraw
    D::Array{Int64,2}
    n::Int64

    # Constructors
    SSAraw(model::SSAmodel) = new(zeros(Int64, model.n_spec, model.t_len), 0)
    SSAraw(D::Array{Int64,2}, n::Int64) = new(D, n)
end

# ? --- SSA Raw Data Functions ---
+(raw1::SSAraw, raw2::SSAraw) = SSAraw(raw1.D .+ raw2.D, raw1.n + raw2.n)

@inline function record!(raw::SSAraw, D, ti)
    @inbounds @simd for k in 1:length(D)
        raw.D[k, ti] += D[k]
    end
end

# ! ----- SSA Actual Data Type -----
struct SSAdata
    D::Array{Float64,2}
    n::Int64

    # Constructors
    SSAdata(raw::SSAraw) = new(raw.D' ./ raw.n, raw.n)
    SSAdata(D::Array{Float64,2}, n::Int64) = new(D, n)
end

Base.zero(::Type{SSAdata}) = SSAdata(zeros(Float64, 0, 0), 0)
Base.iszero(data::SSAdata) = isempty(data.D)

# ? --- SSA Actual Data functions ---
function +(data1::SSAdata, data2::SSAdata)
    if iszero(data1)
        return data2
    elseif iszero(data2)
        return data1
    else
        n1 = data1.n
        n2 = data2.n
        n3 = n1 + n2
        D = similar(data1.D)
        @. D = (data1.D * n1 + data2.D * n2) / n3
        return SSAdata(D, n3)
    end
end

function +(data::SSAdata, raw::SSAraw)
    if iszero(data)
        return SSAdata(raw.D' ./ raw.n, raw.n)
    else
        n1 = data.n
        n2 = raw.n
        n3 = n1 + n2
        D = similar(data.D)
        @. D = (data.D * n1 + raw.D') / n3
        return SSAdata(D, n3)
    end
end