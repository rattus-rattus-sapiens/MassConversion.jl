module MassConversion
import Base.+
using StaticArrays
using Dates
using UUIDs
using JLD2

include("data_structures.jl")
export MCMmodel, MCMoutput

include("mcm_functions.jl")

include("data_io.jl")

end