using Revise
using MassConversion
using Plots

tf = 10;
dt = 1e-2;
IC = [0, 100];
λ = [1e0, 1.0];
# global const tₛ = [(20, 35), (50, 65)]

function R!(S, ri)
    if ri == 1
        S[1] -= 1
    elseif ri == 2
        if S[2] < 1
            r = rand()
            if S[2] > r 
                S[1] += 1
                S[2] = 0 
            end
        else
            S[1] += 1
            S[2] -= 1
        end
    elseif ri == 3
        S[1] -= 1
        S[2] += 1
    else
        throw("invalid reaction attempted - likely negative propensities afoot")
    end
    return nothing
end;

function A!(a, S, t, L)
    a[1] = L[1]*S[1]
    a[2] = L[2]*S[2]*(S[1] + S[2] <= 50)
    a[3] = L[2]*S[1]*(S[1] + S[2] > 50)
end;

function F!(dxdt, S, t, L)
    dxdt[2] = -L[1]*S[2]
end;

rn = 1000;

rec = run_mcm(tf, dt, IC, λ, 3, R!, F!, A!, rn)

plot(rec.mean_discrete)
plot!(rec.mean_continuum)
plot!(rec.mean_total)