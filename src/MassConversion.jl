module MassConversion
import Base.in
import Base.union
import Base.iterate
using Serialize
using Dates
using Plots

include("data_structures.jl")
include("data_functions.jl")
include("mcm_functions.jl")

export run_mcm, run_ssa

end