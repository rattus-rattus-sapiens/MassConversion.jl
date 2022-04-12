function plot(M::MCMdata; discrete::Bool=true, continuous::Bool=true, total::Bool=true)
    p = plot();
    return plot(p, M; discrete=discrete, continuous=continuous, total=total)
end

function plot(S::SSAdata; discrete::Bool=true)
    p = plot();
    return plot(p, S; discrete=discrete)
end

function plot(p::Plots.Plot{B}, M::MCMdata; discrete::Bool=true, continuous::Bool=true, total::Bool=true) where {B<:AbstractBackend}
    if discrete
        plot!(p, M.t_range, M.D; label="Discrete");
    end
    if continuous
        plot!(p, M.t_range, M.C; label="Continuous");
    end
    if total
        plot!(p, M.t_range, M.D .+ M.C; label="Total");
    end
    return p
end

function plot(p::Plots.Plot{B}, S::SSAdata; discrete::Bool=true) where {B<:AbstractBackend}
    if discrete
        plot!(p, S.t_range, S.D; label="SSA");
    end
    return p
end