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

# ! --- MCM Output Type ---
mutable struct MCMraw
    D::Array{Int64,2}
    C::Array{Float64,2}
    n::Int64
    function MCMraw(M::MCMmodel)
        D = zeros(Float64, M.n_spec, M.t_len)
        C = zeros(Float64, M.n_spec, M.t_len)
        new(D, C, 0)
    end
    function MCMraw(D::Array{Int64,2}, C::Array{Float64,2}, n::Int64)
        new(D, C, n)
    end
end

function +(model1::MCMraw, model2::MCMraw)
    MCMraw(model1.D .+ model2.D, model1.C .+ model2.C, model1.n + model2.n)
end

@inline function record!(out::MCMraw, D, C, ti)
    @inbounds @simd for k in 1:length(D)
        out.D[k, ti] += D[k]
        out.C[k, ti] += C[k]
    end
end

# ! --- MCM Data Type ---
struct MCMdata
    D::Array{Float64,2}
    C::Array{Float64,2}
    n::Int64
    function MCMdata(M::MCMraw)
        new(M.D' ./ M.n, M.C' ./ M.n, M.n)
    end
    function MCMdata(D, C, n)
        new(D, C, n)
    end
end

Base.zero(::Type{MCMdata}) = MCMdata(zeros(Float64, 0, 0), zeros(Float64, 0, 0), 0)
Base.iszero(data::MCMdata) = isempty(data.D) || isempty(data.C)

function +(data1::MCMdata, data2::MCMdata)
    if iszero(data1)
        return data2
    elseif iszero(data2)
        return data1
    else
        n1 = data1.n
        n2 = data2.n
        n3 = data1.n + data2.n
        @. D = (data1.D * n1 + data2.D * n2) / n3
        @. C = (data1.C * n1 + data2.C * n2) / n3
        return MCMdata(D, C, n3)
    end
end

function +(data::MCMdata, raw::MCMraw)
    if iszero(data)
        return MCMdata(raw.D' ./ raw.n, raw.C' ./ raw.n, raw.n)
    else
        n1 = data.n
        n2 = raw.n
        n3 = n1 + n2
        @. D = (data.D * n1 + raw.D') / n3
        @. C = (data.C * n1 + raw.C') / n3
        return MCMdata(D, C, n3)
    end
end