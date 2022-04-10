function load_data(path::String)
    files = "$path/" .* readdir(path)
    data = load(files[1], "data") 
    for filename in files[2:end]
        data += load(filename, "data")
    end
    return data
end

function load_data(path::String, samples::Int64)

end