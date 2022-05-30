using LinearAlgebra, Statistics
using BenchmarkTools, Interpolations, LaTeXStrings, Parameters, Plots, QuantEcon, Roots
using Optim, Random

function K!(Kg, g, grid, β, ∂u∂c, f, f', shocks)
    g_func = LinearInterpolation(grid, g, extrapolation_bc = Line())
    for (i, y) in enumerate(grid)
        function h(c)
            vals = ∂u∂c.(g_func.(f(y-c) * shocks)) .* f'(y-c) .* shocks 
            return ∂u∂c - β * mean(vals)
        end
        Kg[i] = find_zero(h, (1e-10, y - 1e-10))
    end
    return Kg
end

K(g, grid, β, ∂u∂c, f, f', shocks) = K!(similar(g), g, grid, β, ∂u∂c, f, f', shocks)

using Optim

function T(w, grid, β, )
    w_func = LinearInterpolation(grid, w, extrapolation_bc = Line())

    if compute_policy
        σ = similar(w)
    end

    for (i, y) in enumerate(grid) 
        objective(c) = u(c) + β * mean(w_func.(f(y-c) .* shocks))
        res = maximize(objective, 1e-10, y)

        if compute_policy
            σ[i] = Optim.maximizer(res)
        end

        Tw[i] = Optim.maximizer(res) 
    end

    if compute_policy
        return Tw, σ
        else
            return Tw
    end
end

isoelastic(c, γ) = isone(γ) ? log(c) : (c^(1 - γ) - 1) / (1 - γ)
Model = @with_kw (α = 0.65,                            # Productivity parameter
                  β = 0.95,                            # Discount factor
                  γ = 1.0,                             # Risk aversion
                  μ = 0.0,                             # First parameter in lognorm(μ, σ)
                  s = 0.1,                             # Second parameter in lognorm(μ, σ)
                  grid = range(1e-6, 4, length = 200), # Grid
                  grid_min = 1e-6,                     # Smallest grid point
                  grid_max = 4.0,                      # Largest grid point
                  grid_size = 200,                     # Number of grid points
                  u = (c, γ = γ) -> isoelastic(c, γ),  # utility function
                  ∂u∂c = c -> c^(-γ),                  # u′
                  f = k -> k^α,                        # production function
                  f′ = k -> α * k^(α - 1),             # f′
)


m = Model();

using Random
Random.seed!(42)

shock_size = 250
shocks = collect(exp.(m.μ .+ m.s * randn(shock_size)));





