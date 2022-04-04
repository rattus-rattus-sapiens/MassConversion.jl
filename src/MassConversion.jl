module MassConversion
import Base.in
import Base.union
import Base.iterate
using JLD2
using Dates
using Plots
using Parameters
using Statistics

include("data_structures.jl")
export to_dict

include("math_functions.jl")

include("mcm_functions.jl")
export run_mcm

include("gillespie_functions.jl")
export run_ssa

include("data_io.jl")
export save_data, load_data

include("plotting.jl")
export plot_total, plot_both, plot_init, plot_save, plot_rel_err, plot_abs_err

end