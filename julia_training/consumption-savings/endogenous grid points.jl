using LinearAlgebra, Statistics
using BenchmarkTools, Interpolations, LaTeXStrings, Parameters, Plots, QuantEcon, Random, Roots

function coleman_egm(g, k_grid, β, u′, u′_inv, f, f′, shocks)

    c = similar(k_grid)

    for (i,k) in enumerate(k_grid)
        vals = u′.(g.(f(k)*shocks)).*f′(k).*shocks
        c[i] = u′_inv(β*mean(vals))
    end

    y = k_grid + c

    Kg = LinearInterpolation(y, c, extrapolation_bc = Line())
    Kg_f(x) = Kg(x)
    return Kg_f
end

function K!(Kg, g, grid, β, u′, f, f′, shocks)
    g_func = LinearInterpolation(grid, g, extrapolation_bc = Line())
    for (i, y) in enumerate(grid)
        function h(c)
            vals = u′.(g_func.(f(y-c) * shocks)) .* f′(y-c) .* shocks 
            return u′(c) - β * mean(vals)
        end
        Kg[i] = find_zero(h, (1e-10, y - 1e-10))
    end
    return Kg
end

K(g, grid, β, u′, f, f′, shocks) = K!(similar(g), g, grid, β, u′, f, f′, shocks)


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
shocks = collect(exp.(mlog.μ .+ mlog.s * randn(shock_size)));

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
    c_star_new = coleman_egm(c_star, k_grid, m.β, m.u′, m.u′, m.f, m.f′, shocks)

    plt = plot()
    plot!(plt, k_grid, c_star.(k_grid), lw = 2, label = L"optimal policy $c^*$")
    plot!(plt, k_grid, c_star_new.(k_grid), lw = 2, label = L"Kc^*")
    plot!(plt, legend = :topleft)
    
end

verify_true_policy(mlog, shocks, c_star)

n = 15
function convergence_check(m, shocks, c_star, g_init, n_iter)
    k_grid = m.grid
    g = g_init
    plt = plot()
    plot!(plt, m.grid, g.(m.grid),
    color = RGBA(0,0,0,1), lw = 2, alpha = 0.6, label = L"initial condition $c(y) = y$")
    for i in 1:n_iter
        new_g = coleman_egm(g, k_grid, m.β, m.u′, m.u′, m.f, m.f′, shocks)
        g = new_g
        plot!(plt, k_grid, new_g.(k_grid), alpha = 0.6, color = RGBA(0,0,(i/n_iter),1), lw = 2, label = "")
    end
    plot!(plt, k_grid, c_star.(k_grid),
    color = :black, lw = 2, alpha = 0.8, label = L"true policy function $c^*$")
    plot!(plt, legend = :topleft)
end

convergence_check(mlog, shocks, c_star, identity, n)


mcrra = Model(α = 0.65, β = 0.95, γ = 1.5)
u′_inv(c) = c^(-1 / mcrra.γ)

crra_coleman(g, m, shocks) = K(g, m.grid, m.β, m.u′, m.f, m.f′, shocks)

crra_coleman_egm(g, m, shocks) = coleman_egm(g, m.grid, m.β, m.u′, u′_inv, m.f, m.f′, shocks)


function coleman(m = m, shocks = shocks; sim_length = 20)
    g = m.grid
    for i in 1:sim_length
        g = crra_coleman(g, m, shocks)
    end
    return g
end

function egm(m, g = identity, shocks = shocks; sim_length = 20)
    for i in 1:sim_length
        g = crra_coleman_egm(g, m, shocks)
    end
    return g.(m.grid)
end

@benchmark coleman($mcrra)
@benchmark egm($mcrra)
