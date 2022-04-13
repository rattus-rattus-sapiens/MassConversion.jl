function load_raw(dataloc::String)
    path = joinpath(pwd(), dataloc)
    dircontents = readdir(path)
    mcm_data = nothing
    ssa_data = nothing
    if "MCM" ∈ dircontents
        files = readdir("$path/MCM", join=true)
        mcm_data = Tuple(MCMdata(load(file, "data_mcm")) for file in files)
    end
    if "SSA" ∈ dircontents
        files = readdir("$path/SSA", join=true)
        ssa_data = Tuple(SSAdata(load(file, "data_ssa")) for file in files)
    end
    println("Loaded")
    return mcm_data, ssa_data
end

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--repno", "-r"
        help = "Number of repeats to simulate"
        arg_type = Int
        "--blocksize", "-b"
        help = "Number of repeats per save datafile"
        arg_type = Int
        "--dir_name", "-d"
        help = "Name of saved data directory"
        arg_type = String
        default = Dates.format(now(), "e-dd-HH:MM:SS")
    end
    return parse_args(s)
end