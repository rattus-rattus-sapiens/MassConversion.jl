using MassConversion
using Plots

const tf = 1.0
const dt = 1e-2
const IC = [1.0,0,0,0]
# λ = [δ, K, k₀, k₁, a₀, a₁]
const λ = [1, 1, 0.5, 0.5, 0, 1e1]
const R = [
    0  1 -1  0
   -1  0  0  1
]

function A!(alpha, state, λ)
    alpha[1] = λ[1]*state[2]
    alpha[2] = λ[2] * λ[3] * (1-state[1])
    alpha[3] = λ[2] * λ[4] * state[1] 
    alpha[4] = (state[1]==0)*λ[5] + (state[1]==1)*λ[6]
    return nothing
end

function F!(dxdt, state, λ)
    dxdt[1] = 0
    dxdt[2] = 0
    return nothing
end

const rn = 100

rec = run_mcm(tf, dt, IC, λ, R, F!, A!, rn)