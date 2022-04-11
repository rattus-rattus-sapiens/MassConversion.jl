module MassConversion
import Base.+, Base.adjoint
using ArgParse
using StaticArrays
using Dates
using UUIDs
using JLD2

abstract type AbstractModel end

include("structures_mcm.jl")
export MCMmodel

include("structures_ssa.jl")
export SSAmodel

include("simulate.jl")
export par_run_sim

include("data_io.jl")
export load_raw, load_ensemble

end