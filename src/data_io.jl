function save_data(rec, casename::String)
    filename = "examples/" * casename * string(floor(now(), Dates.Second)) * ".jld2"
    jldsave(filename; rec)
    println("Saved at " * filename)
end

function save_data(rec, casename::String, filename::String)
    filename = "examples/" * casename  * "/dat/" * filename * ".jld2"
    jldsave(filename; rec)
    println("Saved at " * filename)
end

function load_data(casename::String, filename::String)
    if filename[end-4 : end] == ".jld2"
        jldopen("examples/" * casename * "/dat/" * filename) do file
            return Dict(key => file[key] for key in keys(file)) 
        end
    else
        jldopen("examples/" * casename * "/dat/" * filename * ".jld2", "r") do file
            return Dict(key => file[key] for key in keys(file))
        end
    end
end

function load_data(casename::String, filename::String, keyname::String)
    jldopen("examples/" * casename * "/dat/" * filename) do file
        return file[keyname]
    end
end