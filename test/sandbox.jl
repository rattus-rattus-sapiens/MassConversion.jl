using Catalyst
using MomentClosure
using Latexify
using BenchmarkTools
using DifferentialEquations
import SymbolicUtils.term

rn1 = @reaction_network begin
    (k), X + X → 0
end k

rn2 = @reaction_network begin
    (k), X + X → 0
    (k), X + Y → 0
    (k), Y + X → 0
    (k), Y + Y → 0
end k

eqs1 = generate_raw_moment_eqs(rn1, 2, combinatoric_ratelaw=false)
latexify(eqs1)

eqs2 = generate_raw_moment_eqs(rn2, 2, combinatoric_ratelaw=false)
latexify(eqs2)

rn = @reaction_network begin
    (k_1, k_2), C ↔ 2*C
    (k_1, k_2), D ↔ 2D
    k_2, C + D → C
    k_2, D + C → D
end k_1 k_2 

pmap = (:C => 10, :D => 10)
umap = (:k_1 => 1e1, :k_2 => 1e-2)
tspan = (0., 1.)

dprob = DiscreteProblem(rn, umap, tspan, pmap)
jprob = JumpProblem(rn, dprob, Direct(), save_positions=(false,false))

sol = solve(jprob, SSAStepper(), saveat=10.)
plot(sol)

# !!! here is where the magic happens
function str_to_var(str::String)
    eval(Meta.parse("@variables $str"*"(t)"))
end

function split_var(T::Num)
    return str_to_var(string(T.val.f)*"_c")[1], str_to_var(string(T.val.f)*"_d")[1]
end

function split_var(T::Term)
    return str_to_var(string(T.f)*"_c")[1], str_to_var(string(T.f)*"_d")[1]
end

function split_reac(r::Reaction)
    rate = r.rate
    if sum(r.substoich) == 1
        c, d = split_var(r.substrates[1])
        ps_c = [split_var(p)[1] for p in r.products]
        ps_d = [split_var(p)[2] for p in r.products]
        reacs = [Reaction(rate, [c], ps_c, [1], r.prodstoich), Reaction(rate, [d], ps_d, [1], r.prodstoich)]
        species = Set([d, c, ps_c..., ps_d...])
    end
    return rate, species, reacs
end



# ? ----------------------------------------------- For use later, probably
closed = convert(ODESystem, rn)
closed = generate_raw_moment_eqs(rn, 1, combinatoric_ratelaw=false)
closed = moment_closure(closed, "zero")
latexify(closed)
p = [1e2, 1e3, 1e-2, 9900, 1e3, 1e2]
tspan = (0., 1.)
u = [10, 10, 10]
umap = deterministic_IC(u, closed)
prob = ODEProblem(closed, umap, tspan, p; dt=0.1)
integrator = init(prob, Euler())
h(x) = integrator.f(x, p, 2)

function get_winner()
    r1 = rand(1:6) + rand(1:6)
    r2 = rand(1:6) + rand(1:6)
    while true
        if r1 == 7 & r2 == 7
            return 0
        elseif (r1 == 8) & (r2 == 7)
            return 1
        end
        r1 = r2
        r2 = rand(1:6) + rand(1:6)
    end
end

vec = zeros(100_000_000)

for i = 1:100_000_000
    vec[i] = get_winner()
end 

