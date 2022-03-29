module MassConversion
import Base.in
import Base.union
import Base.iterate

include("data_structures.jl")
include("mcm_functions.jl")

export run_mcm, run_ssa, Interval, UInterval

end