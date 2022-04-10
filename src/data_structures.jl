abstract type AbstractModel end

# TODO: Implement interface with DiffEqs.jl
struct SSAModel <: AbstractModel
end

struct ODEModel <: AbstractModel
end

const TimeRange = 
const Thresholds = Base.Iterators.Enumerate{Vector{Tuple{Float64, Float64}}}

@with_kw struct MCMmodel{F,G} <: AbstractModel
    t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    D_init::Vector{Int64}
    C_init::Vector{Float64}
    λ_reac::Vector{Float64}
    R_mat::Matrix{Int64}
    θ_list::Vector{Tuple{Float64, Float64}}
    rep_num::Int64
    alpha::F
    dxdt::G
end