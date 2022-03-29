using Revise
using MassConversion
using Plots

const tf = 90
const dt = 1e-2
const IC = [4, 0]
const λ = [1e-3, 6e-1, 0, 1]
global const tₛ = [(20, 35), (50, 65)]

const R = [
   -1  0  0  1  1 -1
    0 -1 -1  0 -1  1
]

function A!(a, S, t, L)
    a[1] = L[1]*S[1]*(S[1]-1)
    a[2] = L[1]*S[1]*S[2]
    a[3] = L[1]*S[1]*S[2]
    for pair in tₛ
        if pair[1] < t < pair[2]
            a[4] = L[3]*S[1]
            break
        else
            a[4] = L[2]*S[1]
        end
    end
    a[5] = L[4]*S[2]*(S[1] + S[2] <= 150)
    a[6] = L[4]*S[1]*(S[1] + S[2] > 150)
end

function F!(dxdt, S, t, L)
    for pair in tₛ
        if pair[1] < t < pair[2]
            dxdt[2] = S[2]*(L[3] - L[1]*S[2])
            break
        else
            dxdt[2] = S[2]*(L[2] - L[1]*S[2])
        end
    end
end

const rn = 1

rec = run_mcm(tf, dt, IC, λ, R, F!, A!, rn)

tspan = 0:dt:tf

plot(tspan, rec.mean_discrete)
plot!(tspan, rec.mean_continuum)
plot!(tspan, rec.mean_total)