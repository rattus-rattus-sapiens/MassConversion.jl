using MassConversion
using Base.Threads
using BenchmarkTools
using Serialization
using Dates

const rn = 100# Meta.parse(ARGS[1])

function main(is_ssa::Bool, rn::Int64)
    tf = 80
    dt = 1e-3
    λ = [1e-3, 6e-1, 1]
    R = [
        
    ]
end