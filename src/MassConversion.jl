module MassConversion
import Base.in
import Base.union
import Base.iterate
using JLD2
using Dates
using Plots

include("data_structures.jl")
include("data_manipulation.jl")
include("mcm_functions.jl")
include("gillespie_functions.jl")

export run_mcm, run_ssa

end