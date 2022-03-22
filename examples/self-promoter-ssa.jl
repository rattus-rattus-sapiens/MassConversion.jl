using Revise
using MassConversion
using Plots

const tf = 20
const dt = 1e-1
const IC = [1000,0,0]
# λ = [Λ, θ, K, β, δ, α₀, α₁]
const λ = [0.15, 10000, 0.15, 10, 1, 50, 1000]
const R = [
   -2  2  0  0 -1  1
    1 -1 -1  1  0  0
    0  0  1 -1  0  0
]

function A!(alpha, state, λ)
    alpha[1] = λ[1]*state[1]*(state[1]-1)
    alpha[2] = λ[1]*λ[2]*state[2]
    alpha[3] = λ[3]*state[2]*(1-state[3])
    alpha[4] = λ[3]*λ[4]*state[3]
    alpha[5] = λ[5]*state[1]
    alpha[6] = λ[6]*(1-state[3]) + λ[7]*state[3]
    return nothing
end

const rn = 1

rec = run_ssa(tf, dt, IC, λ, R, A!, rn)

Plots.plot(rec.mean_total')