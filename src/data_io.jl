function save_data(data::Any)
    filename = string(floor(now(), Dates.second)) * ".jld2"
    jldsave(filename, data)
end

function save_data(data, casename::String)
    filename = "examples/" * casename * string(floor(now(), Dates.Second)) * ".jld2"
    jldsave(filename; data)
    println("Saved at " * filename)
end

function save_data(data, casename::String, filename::String)
    filename = "examples/" * casename  * "/dat/" * filename * ".jld2"
    jldsave(filename; data)
    println("Saved at " * filename)
end

function load_data(casename::String, filename::String)
    if filename[end-4 : end] == ".jld2"
        jldopen("examples/" * casename * "/dat/" * filename) do file
            return file["data"]
        end
    else
        jldopen("examples/" * casename * "/dat/" * filename * ".jld2", "r") do file
            return file["data"]
        end
    end
end

function load_data(casename::String, filename::String, keyname::String)
    jldopen("examples/" * casename * "/dat/" * filename) do file
        return file[keyname]
    end
end

function rand_ID()
    str = ""
    for i=1:16
        str*=rand("1234567890ABCDEF")
    end
    return str
end