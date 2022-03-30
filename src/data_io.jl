function save_mcm(rec::O, str::String="test_data") where O
    filename = "./examples/dat/" * str * "-" * string(floor(now(), Dates.Second)) * ".jld2"
    jldsave(filename; rec)
    println("Saved at " * filename)
end

function load_mcm(str::String)
    jldopen(str, "r") do file
        return file["rec"]
    end
end