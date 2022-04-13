module MassConversion
using ArgParse
using StaticArrays
using Dates
using UUIDs
using JLD2

# Base overloads
import Base.+, Base.adjoint, Base.length

# Plots overloads
using Plots
import Plots.plot

abstract type AbstractModel end
abstract type AbstractData end
abstract type AbstractError end
include("structures_mcm.jl")
include("structures_ssa.jl")
include("structures_err.jl")
export MCMmodel, SSAmodel, MCMdata, SSAdata, RelativeError, rep_no

include("simulate.jl")
export par_run_sim

include("data_io.jl")
export load_raw, parse_cmd

include("plotting.jl")
export plot

end