using LinearAlgebra, Statistics
using BenchmarkTools, Interpolations, LaTeXStrings, Parameters, Plots, QuantEcon, Random, Roots

function colemand_egm(g, k_grid, β, u', u'_inv, f, f', shocks)

    c = similar(k_grid)

    for (i,k) in enumerate(k_grid)
        vals = u'.(g.(f(k)*shocks)).*f'(k).*shocks
        c[i] = u'_inv(β*mean(vals))
    end

    y = k_grid + c

    Kg = LinearInterpolation(y, c, extrapolation_bc = Line())
    Kg_f(x) = Kg(x)
    return Kg_f
end

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


Model = @with_kw (α = 0.65, # productivity parameter
                  β = 0.95, # discount factor
                  γ = 1.0,  # risk aversion
                  μ = 0.0,  # lognorm(μ, σ)
                  s = 0.1,  # lognorm(μ, σ)
                  grid_min = 1e-6, # smallest grid point
                  grid_max = 4.0,  # largest grid point
                  grid_size = 200, # grid size
                  u = γ == 1 ? log : c->(c^(1-γ)-1)/(1-γ), # utility function
                  u′ = c-> c^(-γ), # u'
                  f = k-> k^α, # production function
                  f′ = k -> α*k^(α-1), # f'
                  grid = range(grid_min, grid_max, length = grid_size) # grid
)

mlog = Model();

Random.seed!(42)
shock_size = 250 
shocks = exp.(mlog.μ + mlog.s * randn(shock_size));

c_star(y) = (1 - mlog.α * mlog.β) * y

# some useful constants
ab = mlog.α * mlog.β
c1 = log(1 - ab) / (1 - mlog.β)
c2 = (mlog.μ + mlog.α * log(ab)) / (1 - mlog.α)
c3 = 1 / (1 - mlog.β)
c4 = 1 / (1 - ab)

v_star(y) = c1 + c2 * (c3 - c4) + c4 * log(y)

function verify_true_policy(m, shocks, c_star)
    k_grid = m.grid
    c_star_new = colemand_egm(c_star, k_grid, m.β, m.u′, m.u′, m.f, m.f′, shocks)

    plt = plot()
    plot!(plt, k_grid, c_star.(k_grid), lw = 2, label = L"optimal policy $c^*$")
    plot!(plt, k_grid, c_star_new.(k_grid), lw = 2, label = L"Kc^*")
    plot!(plt, legend = :topleft)
    
end
