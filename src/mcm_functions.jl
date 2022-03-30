"""
    sample_dist(a::Union{AbstractArray{T}, Tuple{Vararg{T}}}, a₀::T) where T <: Real

Samples from an unnormalised discrete distribution.

### Arguments
- `a`: Array/Tuple of weights
- `a₀`: Sum of weights
"""
function sample_dist(a::Union{AbstractArray{T},Tuple{Vararg{T}}}, a₀::T) where {T<:Real}

    i = 1
    cumsum = a[1]
    r = rand()

    while cumsum < r * a₀
        i += 1
        cumsum += a[i]
    end

    return i
end

function exec_transition!(state::Vector, Kn::Int64, tid::Int64)
    if tid <= Kn
        state[tid] -= 1
        state[tid+Kn] += 1
    else
        if state[tid] < 1
            if rand() < state[tid]
                state[tid] = 0
                state[tid-Kn] += 1
            end
        else
            state[tid] -= 1
            state[tid-Kn] += 1
        end
    end
end

function B!(_alpha, _state, λ, enumθ, Rn, Kn)
    for (j, pair) in enumθ
        _alpha[Rn+j] = λ[Rn+j] * _state[j] * (_state[j] + _state[j+Kn] > pair[2])
        _alpha[Rn+j+Kn] = λ[Rn+j] * _state[j+Kn] * (_state[j] + _state[j+Kn] < pair[1])
    end
    return nothing
end

"""
    run_mcm(tf, dt, IC, λ, R, F!, A!, rn)

Run the mass-conversion method.

### Arguments
- `tf` : Simulation run time
- `dt` : Time step
- `IC` : Vector containing initial condition
- `λ` : Vector containing rate parameters. The last K parameters, where K is the number of species, is assumed to be the regime transition rate.
- `R` : Stoichiometric matrix representing all non-transitional reactions
- `θ` : Vector containing tuples representing the lower and upper thresholds, respectively
- `F!` : Function for computing dxdt for the forward Euler method
- `A!` : Function for computing propensity functions for all non-transitional reactions
- `rn` : Number of repeats
"""
function run_mcm(
    tf, 
    dt, 
    IC::Vector{T}, 
    λ::Vector{Float64}, 
    R::Matrix{Int64}, 
    θ::Vector{Tuple{S,P}}, 
    F!::Function, 
    A!::Function, 
    rn::Int64
    ) where {T<:Real,S<:Real,P<:Real}

    # Get constants
    tn = floor(Int64, tf / dt) + 1
    Kn = floor(Int64, length(IC) / 2) # Number of unique species
    Kt = 2 * Kn
    Rn = size(R, 2) # Number of non-transitional reactions

    # Input validation
    if length(IC) / 2 ≠ Kn
        throw("even number of input species expected")
    end
    if size(R, 1) ≠ length(IC)
        throw("ill-defined stoichimetric matrix")
    end
    if length(θ) ≠ Kn
        throw("number of specified thresholds mismatched with number of species")
    end

    # Input sanitisation
    IC = Vector{Float64}(IC)
    θ = [Tuple(Float64[θ[i]...]) for i in 1:length(θ)] # TODO?: Casts tuple type to (Float64, Float64)

    # Preallocation
    rec = zeros(2 * Kn, tn, rn)
    state = similar(IC)
    α = zeros(Rn + 2 * Kn)
    dxdt = zeros(2 * Kn)

    # Propensity function for transitional reactions
    enumθ = enumerate(θ)

    for ri ∈ 1:rn
        # Initial conditions
        t = 0.0
        td = dt
        state .= IC
        ti = 1
        rec[:, ti, ri] .= state

        while t < tf
            # Update propensity functions
            A!(α, state, t, λ)
            B!(α, state, λ, enumθ, Rn, Kn)
            α₀ = sum(α)

            # Time to next reaction
            τ = log(1 / rand()) / α₀

            if t + τ < td
                # Execute stochastic event
                t += τ
                reaci = sample_dist(α, α₀)
                if reaci <= Rn
                    for i in 1:Kt
                        state[i] += R[i, reaci]
                    end
                else
                    exec_transition!(state, Kn, reaci - Rn)
                end
            else
                # Execute ODE update
                F!(dxdt, state, t, λ)
                @. state += dt * dxdt

                # Record
                ti += 1
                rec[:, ti, ri] .= state
                t = (ti - 1) * dt
                td = ti * dt
            end
        end

        # Final record
        rec[:, end, ri] .= state
    end

    return MCMOutput(rec, tf, dt, IC, λ, R, θ, rn)
end

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