module MassConversion
import Base.+, Base.adjoint
using ArgParse
using StaticArrays
using Dates
using UUIDs
using JLD2

# overloading plot here - don't worry - no type piracy!
using Plots
import Plots.plot

abstract type AbstractModel end

include("structures_mcm.jl")
export MCMmodel

include("structures_ssa.jl")
export SSAmodel

include("simulate.jl")
export par_run_sim

include("data_io.jl")
export load_raw, load_ensemble

include("plotting.jl")
export plot

end