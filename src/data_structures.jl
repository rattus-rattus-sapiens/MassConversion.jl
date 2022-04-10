abstract type AbstractModel end

# TODO: Implement interface with some kind of Gillespie package (DiffEqs possible?)
struct SSAModel <: AbstractModel end

# TODO: Implement interface with DiffEqs.jl
struct ODEModel <: AbstractModel end

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
        calc_dxdt::Function) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

        # get derived params
        t_start = Float64(t_span.ref)
        t_step = Float64(t_span.step)
        t_len = t_span.len
        t_final = t_start + t_step * (t_len-1)
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
mutable struct MCMoutput
    D::Array{Int64,2}
    C::Array{Float64,2}
    n::Int64
    function MCMoutput(M::MCMmodel)
        D = zeros(Float64, M.n_spec, M.t_len)
        C = zeros(Float64, M.n_spec, M.t_len)
        new(D, C, 0)
    end
    function MCMoutput(D::Array{Int64, 2}, C::Array{Float64, 2}, n::Int64)
        new(D,C,n)
    end
end

function +(model1::MCMoutput, model2::MCMoutput)
    MCMoutput(model1.D .+ model2.D, model1.C .+ model2.C, model1.n + model2.n) 
end

@inline function record!(out::MCMoutput, D, C, ti)
    @inbounds @simd for k in 1:length(D)
        out.D[k, ti] += D[k]
        out.C[k, ti] += C[k]
    end
end