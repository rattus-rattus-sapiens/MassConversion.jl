function load_raw(path::String)
    files = "$path/" .* readdir(path)
    n_files = length(files)
    println("Loading $n_files files...")
    data = []
    for filename in files
        data += load(filename, "data")
    end
    println("Files loaded")
    return MCMdata(data)
end

function load_ensemble(path::String, samples::Int64)
    files = "$path/" .* readdir(path)
    n_files = length(files)
    files_per_block = n_files ÷ samples
    println("Loading $samples batches of $files_per_block files")
    data = zeros(MCMdata, samples)
    for j = 1:samples
        for filename in files[(j-1)*files_per_block+1 : j*files_per_block]
            data[j] += load(filename, "data")
        end
    end
    return data
end