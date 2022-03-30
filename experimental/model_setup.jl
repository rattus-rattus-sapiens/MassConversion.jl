import Base.<=
    
abstract type AbstractModel end
abstract type AbstractSpecies end
abstract type AbstractReaction end

struct ModelSpecError <: Exception
    val::String
    str::String
end

Base.showerror(io::IO, e::ModelSpecError) = print(io, e.var, e.str)


struct Model <: AbstractModel
    S::Vector{String}
    D::OrderedDict{String, Int64}
    R::Matrix{Int64}
    function Model()
        S = Vector{String}()
        D = OrderedDict{String, Int64}()
        R = Matrix{Int64}(undef, (0, 0))
        new(S,D,R)
    end
end

function get_init(m::AbstractModel) m.D end

struct Species <: AbstractSpecies
    s::String
    init::Int64
    function Species(s::String, init::Int64)
        try
            s::Symbol = Meta.parse(s)
        catch
            throw(ModelSpecError(s, " is an invalid species identifier"))
        end
        new(string(s), init)
    end
end

struct Reaction <: AbstractReaction
    k::Float64
    reactants::OrderedDict{String, Int64}
    products::OrderedDict{String, Int64}
    function Reaction(k::Float64, s::String)
        reactants, products = reaction_parse(s)
        new(k, reactants, products)
    end
end

function <=(m::AbstractModel, s::AbstractSpecies)
    if s.s ∈ m.S
        throw(ModelSpecError(s.s, " already initialised"))
    else
        m.D[s.s] = s.init
        push!(m.S, s.s)
        nothing
    end
end

function reaction_parse(str::String)
    reactants = OrderedDict{String, Int64}()
    products = OrderedDict{String, Int64}()
    expr = Meta.parse(str)
    if expr.head == :-->
        add_species(reactants, expr.args[1])
        add_species(products, expr.args[2]) 
    else
        throw("Ill-posed reaction passed")
    end
    return reactants, products
end

function add_species(d::AbstractDict, expr::Expr)
    if expr.args[1] == :*
        n = get(d, string(expr.args[3]), 0)
        d[string(expr.args[3])] = n + expr.args[2]
    else
        for i in 2:length(expr.args)
            add_species(d, expr.args[i])
        end
    end
end

function add_species(d::AbstractDict, sym::Symbol)
    n = get(d, string(sym), 0)
    d[string(sym)] = n + 1
end