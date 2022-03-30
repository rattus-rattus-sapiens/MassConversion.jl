"""
    record_state!()
"""
function record_state!(rec::Array, t::Float64, tn::Int64, tp::Array{Int64}, Δt::Float64, state::Array, repi::Int64)
    if t > tp[1] * Δt
        _tc::Int64 = min(tn, floor(t / Δt) + 1)
        for i ∈ 1:length(state)
            for j ∈ tp[1]:_tc
                rec[i, j, repi] += state[i]
            end
        end
        tp .= _tc + 1
    end
    nothing
end

"""
    run_ssa(tf, dt, IC, λ, R, A!, rn)
"""
function run_ssa(tf, dt, IC::Vector{T}, λ::Vector{Float64}, R::Matrix, A!::Function, rn::Int64) where {T<:Real}

    # Get timesteps etc
    tn = floor(Int64, tf / dt) + 1
    K = length(IC)

    # Preallocation
    rec = zeros(K, tn, rn)
    state = similar(IC)
    α = zeros(size(R, 2))
    reac_count = zeros(Int64, size(R, 2))

    for ri ∈ 1:rn
        # Initial conditions
        t = 0.0
        tp = [1]
        state .= IC

        while t < tf
            # Update propensity functions
            A!(α, state, λ)
            α₀ = sum(α)

            # Time to next reaction
            τ = log(1 / rand()) / α₀

            # Execute stochastic event
            t += τ
            reaci = sample_dist(α, α₀)
            for i ∈ 1:K
                state[i] += R[i, reaci]
            end
            reac_count[reaci] += 1

            # Record
            record_state!(rec, t, tn, tp, dt, state, ri)
        end

        # Final record
        rec[:, end, ri] .= state
    end

    return SSAOutput(rec)
end