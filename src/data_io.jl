function save_mcm(rec, str::String="test_data")
    filename = "./examples/dat/" * str * "-" * string(floor(now(), Dates.Second)) * ".jld2"
    jldsave(filename; rec)
    println("Saved at " * filename)
end

function load_mcm(str::String)
    jldopen(str, "r") do file
        return file["rec"]
    end
end