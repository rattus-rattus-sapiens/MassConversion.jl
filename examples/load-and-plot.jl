using FileIO, JLD2 # Data loading dependencies
using Plots

# ! This must be set to the correct data location to load data 
loadstr = "examples/dat/alt-log-2022-03-29T01:28:52.jld2"

# * Load data
rec = FileIO.load(loadstr, "rec")