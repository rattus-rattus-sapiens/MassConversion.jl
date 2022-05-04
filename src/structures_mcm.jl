# ! --- MCM Model Type ---
struct MCMmodel{F,G,S1,S2,S3} <: AbstractModel
    t_start::Float64
    t_step::Float64
    t_final::Float64
    t_len::Int64
    n_spec::Int64
    D_init::SVector{S1,Int64}
    C_init::SVector{S1,Float64}
    n_reac::Int64
    R_mat::Matrix{Int64}
    λ_reac::SVector{S2,Float64}
    λ_tran::SVector{S3,Float64}
    θ::Vector{Tuple{Float64,Float64}}
    calc_prop::F
    calc_dxdt::G

    # Constructor function
    function MCMmodel(
        t_span::AbstractRange,
        D_init::AbstractVector{Int64},
        C_init::AbstractVector{T1},
        λ_reac::AbstractVector{T2},
        λ_tran::Union{AbstractVector{T3},T3},
        R_mat::Matrix{Int64},
        θ::Vector{Tuple{T4,T4}},
        calc_prop::Function,
        calc_dxdt::Function
    ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

        # get derived params
        t_start = Float64(t_span.ref)
        t_step = Float64(t_span.step)
        t_len = t_span.len
        t_final = t_start + t_step * (t_len - 1)
        n_spec = length(D_init)
        n_reac = size(R_mat, 2)

        # input validation
        if length(D_init) ≠ length(C_init)
            throw("number of species in D and C must be equal")
        end
        if size(R_mat, 1) ≠ 2 * n_spec
            throw("ill-formed stoichiometry matrix")
        end
        if length(θ) ≠ n_spec
            throw("insufficent thresholds given")
        end

        # typecast remaining inputs
        D_init = SVector{n_spec,Int64}(D_init)
        C_init = SVector{n_spec,Float64}(C_init)
        λ_reac = SVector{n_reac,Float64}(λ_reac)

        if λ_tran isa Real
            λ_tran = SVector{2 * n_spec,Float64}(λ_tran .* ones(2 * n_spec))
        else
            λ_tran = SVector{2 * n_spec,Float64}(λ_tran)
        end

        new{typeof(calc_prop),typeof(calc_dxdt),n_spec,n_reac,2 * n_spec}(
            t_start, t_step, t_final, t_len, n_spec, D_init, C_init, n_reac, R_mat, λ_reac, λ_tran, θ, calc_prop, calc_dxdt
        )
    end
end

get_t_range(M::MCMmodel) = M.t_start : M.t_step : M.t_final

# ! --- MCM Raw Data Type ---
"""
    MCMraw{T}

This function is for things
"""
mutable struct MCMraw{T}
    D::Array{Int64,2}
    C::Array{Float64,2}
    t_range::T
    n::Int64
    function MCMraw(M::MCMmodel)
        D = zeros(Float64, M.n_spec, M.t_len)
        C = zeros(Float64, M.n_spec, M.t_len)
        tr = get_t_range(M)
        new{typeof(tr)}(D, C, tr, 0)
    end
end

@inline function record!(out::MCMraw, D, C, ti)
    @inbounds @simd for k in 1:length(D)
        out.D[k, ti] += D[k]
        out.C[k, ti] += C[k]
    end
end

# ! --- MCM Data Type ---
struct MCMdata{T} <: AbstractData
    D::Array{Float64,2}
    C::Array{Float64,2}
    t_range::T
    n::Int64
    function MCMdata(R::MCMraw)
        tr = R.t_range
        new{typeof(tr)}(R.D' ./ R.n, R.C' ./ R.n, tr, R.n)
    end
    MCMdata(D, C, tr, n) = new{typeof(tr)}(D, C, tr, n)
end

# * --- Getters ---
get_discrete(M::MCMdata) = M.D;
get_continuous(M::MCMdata) = M.C;
get_total(M::MCMdata) = M.D .+ M.C;

# * --- Overloaded functions ---
Base.zero(::Type{MCMdata}) = MCMdata(zeros(Float64, 0, 0), zeros(Float64, 0, 0), 0.0:0.1:0.0, 0)
Base.iszero(data::MCMdata) = isempty(data.D) || isempty(data.C)
+(d::MCMdata) = d
function +(d1::MCMdata, d2::MCMdata)
    if iszero(d1)
        return d2
    elseif iszero(d2)
        return d1
    else
        D = similar(d1.D)
        C = similar(d1.C)
        n1 = d1.n
        n2 = d2.n
        n3 = d1.n + d2.n
        @. D = (d1.D * n1 + d2.D * n2) / n3
        @. C = (d1.C * n1 + d2.C * n2) / n3
        tr = d1.t_range
        return MCMdata(D, C, tr, n3)
    end
end