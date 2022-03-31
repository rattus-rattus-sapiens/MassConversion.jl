module MassConversion
import Base.in
import Base.union
import Base.iterate
using JLD2
using Dates
using Plots
using Parameters

include("data_structures.jl")

include("data_manipulation.jl")
export to_dict

include("mcm_functions.jl")
export run_mcm

include("gillespie_functions.jl")
export run_ssa

include("data_io.jl")
export save_mcm, load_mcm

include("plotting.jl")
export plot_total, plot_both, plot_init

end