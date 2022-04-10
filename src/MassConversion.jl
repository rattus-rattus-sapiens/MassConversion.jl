module MassConversion
import Base.in
import Base.union
import Base.iterate
using JLD2
using Dates
using Plots
using Statistics

include("data_structures.jl")

include("mcm_functions.jl")

include("data_io.jl")
export save_data, load_data, scratch_save

include("plotting.jl")
export plot_total, plot_both, plot_init, plot_save, plot_rel_err, plot_abs_err

end