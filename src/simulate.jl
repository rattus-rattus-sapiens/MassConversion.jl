@inline function sample_dist(a::AbstractVector{T}, a₀::T) where {T<:Real}
    i = 1
    cumsum = a[1]
    r = rand()
    @inbounds while cumsum < r * a₀
        i += 1
        cumsum += a[i]
    end
    return i
end

@inline function move!(D, C, direc::Int64, idx::Int64)
    @inbounds if direc == 1
        D[idx] -= 1
        C[idx] += 1
    elseif direc == -1
        if C[idx] < 1
            if rand() < C[idx]
                D[idx] += 1
                C[idx] = 0
            end
        else
            D[idx] += 1
            C[idx] -= 1
        end
    end
    return nothing
end

@inline function calc_tran(
    fwd::AbstractVector{Float64},
    bwd::AbstractVector{Float64},
    D::AbstractVector{Int64},
    C::AbstractVector{Float64},
    θ::AbstractVector{Tuple{Float64,Float64}},
    L::AbstractVector{Float64})
    @inbounds for (i, pair) in enumerate(θ)
        fwd[i] = L[i] * D[i] * (D[i] + C[i] > pair[2])
        bwd[i] = L[i] * C[i] * (D[i] + C[i] < pair[1])
    end
    return nothing
end

@inline function execute_reaction!(D, C, model::MCMmodel, rid)
    if rid <= model.n_reac
        @inbounds @simd for k = 1:model.n_spec
            D[k] += model.R_mat[k, rid]
            C[k] += model.R_mat[model.n_spec+k, rid]
        end
    elseif rid <= model.n_reac + model.n_spec
        move!(D, C, 1, rid - model.n_reac)
    else
        move!(D, C, -1, rid - model.n_reac - model.n_spec)
    end
end

@inline function execute_reaction!(D, model::SSAmodel, rid)
    @inbounds @simd for k = 1:model.n_spec
        D[k] += model.R_mat[k, rid]
    end
end

@inline function _sim_block(model::MCMmodel{F,G,S1,S2,S3}, blocksize::Int64) where {F,G,S1,S2,S3}
    # malloc 
    data = MCMraw(model)
    a = MVector{S2 + S3,Float64}(zeros(Float64, S2 + S3))
    dxdt = MVector{S1,Float64}(zeros(Float64, S1))
    D = MVector{S1,Int64}(zeros(Int64, S1))
    C = MVector{S1,Float64}(zeros(Float64, S1))

    # Create views
    prop_reac = view(a, 1:model.n_reac)
    prop_tfwd = view(a, model.n_reac+1:model.n_spec)
    prop_tbwd = view(a, model.n_reac+model.n_spec+1:model.n_reac+2*model.n_spec)
    calc_prop = model.calc_prop
    calc_dxdt = model.calc_dxdt

    for _ in 1:blocksize
        # Initial conditions
        D .= model.D_init
        C .= model.C_init
        t = 0.0
        td = model.t_step
        ti = 1

        record!(data, D, C, ti)

        while t < model.t_final
            # Update propensities
            calc_prop(prop_reac, D, C, t, model.λ_reac)
            calc_tran(prop_tfwd, prop_tbwd, D, C, model.θ, model.λ_tran)
            a0 = sum(a)

            # Get next reaction time
            @fastmath τ = log(1 / rand()) / a0

            if t + τ < td
                t += τ
                rid = sample_dist(a, a0)
                execute_reaction!(D, C, model, rid)
            else
                calc_dxdt(dxdt, D, C, t, model.t_step, model.λ_reac)

                # record
                ti += 1
                record!(data, D, C, ti)
                t = (ti - 1) * model.t_step
                td = ti * model.t_step
            end

        end
        data.n += 1
    end

    return data

end

@inline function _sim_block(model::SSAmodel{F,S1,S2}, blocksize::Int64) where {F,S1,S2}
    # malloc
    data = SSAraw(model)
    a = MVector{S2,Float64}(zeros(Float64, S2))
    D = MVector{S1,Int64}(zeros(Int64, S1))

    # create views
    calc_prop = model.calc_prop

    for _ in 1:blocksize
        # initial conditions
        D .= model.D_init
        t = model.t_start
        td = model.t_step
        ti = 1

        record!(data, D, ti)

        while t < model.t_final
            # update propensities
            calc_prop(a, D, t, model.λ_reac)
            a0 = sum(a)

            # get next reaction time
            @fastmath τ = log(1 / rand()) / a0

            if t + τ < td
                t += τ
                rid = sample_dist(a, a0)
                execute_reaction!(D, model, rid)
            else
                # record
                ti += 1
                record!(data, D, ti)
                t = (ti - 1) * model.t_step
                td = ti * model.t_step
            end
        end
        data.n += 1
    end

    return data

end

function par_run_sim(
    model::M,
    rep_no::Int64,
    blocksize::Int64;
    dir_name::String=Dates.format(now(), "e-dd-HH:MM:SS")) where {M<:AbstractModel}

    num_mbs = round((8 * 2 * model.n_spec * model.t_len * rep_no / blocksize) / (1024^2), digits=3)
    n_blocks::Int64 = rep_no / blocksize

    # setup saving directory
    path = joinpath("dat", dir_name)
    if num_mbs > 100
        throw("disc usage will exceed 100 MiBs - estimated $num_mbs MiBs")
    else
        println("Estimated disc usage $num_mbs MiBs. Saving to $path")
        println("Creating $n_blocks files each of size $(num_mbs*1000) KiB")
        mkpath(path)
    end

    if M <: MCMmodel
        O = MCMdata
    elseif M <: SSAmodel
        O = SSAdata
    end

    Threads.@threads for n in 1:n_blocks
        t = @elapsed data = _sim_block(model, blocksize)
        jldsave(joinpath(path, string(uuid4()) * ".jld2"); data=O(data))
        t = round(t, digits=5)
        println("Block $n created, elapsed $t seconds")
    end
end