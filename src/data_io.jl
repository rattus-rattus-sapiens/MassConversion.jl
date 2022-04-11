function load_raw(dataloc::String)
    path = joinpath(pwd(), dataloc)
    dircontents = readdir(path)
    mcm_data = nothing
    ssa_data = nothing
    if "MCM" ∈ dircontents
        files = readdir("$path/MCM", join=true)
        mcm_data = [MCMdata(load(file, "data_mcm")) for file in files]
    end
    if "SSA" ∈ dircontents
        files = readdir("$path/SSA", join=true)
        ssa_data = [SSAdata(load(file, "data_ssa")) for file in files]
    end
    println("Loaded")
    return mcm_data, ssa_data
end