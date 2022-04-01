function update_mu(μ, n, x)
    
end

function update_var(μ, σ, n, x)
    μ = (1 / (n + 1)) * (n * μ + x)
    s = ((n + 1) / n) * ((n / (n + 1)) * σ + (1 / n) * (x - μ)^2)
    σ = (n / (n + 1)) * σ + (1 / n) * (x - μ)^2
    return μ, σ, s
end