push!(LOAD_PATH,"..\\src\\")
using Documenter, MassConversion

makedocs(
    sitename="MassConversion.jl",
    modules = [MassConversion],
    pages = ["Index" => "index.md"]
)

deploydocs(
    repo = "github.com/jkynaston/MassConversion.jl"
)