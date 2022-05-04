function load_raw(dataloc::String)
    path = joinpath(pwd(), dataloc)
    dircontents = readdir(path)
    files = readdir(path, join=true)
    data = Tuple(load(file, "data") for file in files)
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