"""
    run_mcm(tf, dt, IC, λ, R, F!, A!, rn)
"""
function run_mcm(tf, dt, IC::Vector{Float64}, λ::Vector{Float64}, R::Matrix, F!::Function, A!::Function, rn::Int64)

    # Get timesteps etc
    tn = floor(Int64, tf / dt) + 1
    K = length(IC)

    # Preallocation
    rec = zeros(K, tn, rn)
    state = similar(IC)
    α = zeros(size(R, 2))
    dxdt = zeros(K)

    for ri ∈ 1:rn
        # Initial conditions
        t = 0.0
        td = dt
        state .= IC
        ti = 1
        rec[:, ti, ri] .= state

        while t < tf
            # Update propensity functions
            A!(α, state, λ)
            α₀ = sum(α)

            # Time to next reaction
            τ = log(1 / rand()) / α₀

            if t + τ < td
                # Execute stochastic event
                t += τ
                reaci = sample_dist(α, α₀)
                for i ∈ 1:K state[i] += R[i, reaci] end
            else
                # Execute ODE update
                F!(dxdt, state, λ)
                @. state += dt * dxdt

                # Record
                t = (ti - 1) * dt
                rec[:, ti, ri] .= state
                ti += 1
                td = ti * dt
            end
        end

        # Final record
        rec[:, end, ri] .= state
    end

    return rec
end