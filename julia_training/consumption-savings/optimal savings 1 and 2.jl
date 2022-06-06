
## solving consumption savings model, using the time iteration algorithm 

using LinearAlgebra, Statistics
using BenchmarkTools, Interpolations, LaTeXStrings, Parameters, Plots, QuantEcon, Roots
using Optim, Random

using BenchmarkTools, Interpolations, Parameters, Plots, QuantEcon, Roots


### Example optimal growth, ramsey
using Pkg
Pkg.add("QuantEcon")
Pkg.add("Interpolations")
Pkg.add("NLsolve")
Pkg.add("Random")
Pkg.add("Parameters")
Pkg.add("LaTeXStrings")
Pkg.add("BenchmarkTools")
Pkg.add("Roots")    
using LinearAlgebra, Statistics
using Plots, QuantEcon, Interpolations, NLsolve, Optim, Random, Parameters, LaTeXStrings
using Optim: maximum, maximizer 

###

f(x) = 2 .* cos.(6x) .+ sin.(14x) .+ 2.5
c_grid = 0:.2:1
f_grid = range(0,  1, length = 150)

Af = LinearInterpolation(c_grid, f(c_grid))

plt = plot(xlim = (0,1), ylim = (0,6))
plot!(plt, f, f_grid, color = :blue, lw = 2, alpha = 0.8, label = "true function")
plot!(plt, f_grid, Af.(f_grid), color = :green, lw = 2, alpha = 0.8,
      label = "linear approximation")
plot!(plt, f, c_grid, seriestype = :sticks, linestyle = :dash, linewidth = 2, alpha = 0.5,
      label = "")
plot!(plt, legend = :top)

### 


# Bellman operator 

function T(w;p, tol = 1e-10)
    @unpack β, u, f, ξ, y = p # unpack parameters
    w_func = LinearInterpolation(y, w)

    Tw = similar(w)
    σ = similar(w)
    for (i, y_val) in enumerate(y)
        # solve maximization for each point in y, using y itself as initial condition.
        results = maximize(c -> u(c;p) + β * mean(w_func.(f(y_val - c;p) .* ξ)), tol, y_val)
        Tw[i] = maximum(results)
        σ[i] = maximizer(results)
    end
    return (;w = Tw, σ) # returns named tuple of results
end

# Example 

Random.seed!(42) # for reproducible results
u(c;p) = log(c) # utility
f(k;p) = k^p.α # deterministic part of production function
OptimalGrowthModel = @with_kw (α = 0.4, β = 0.96, μ = 0.0, s = 0.1,
                  u = u, f = f, # defaults defined above
                  y = range(1e-5, 4.0, length = 200), # grid on y
                  ξ = exp.(μ .+ s * randn(250)) # monte carlo shocks
) # named tuples defaults

# True value and policy function
function v_star(y;p)
    @unpack α, μ, β = p
    c1 = log(1 - α * β) / (1 - β)
    c2 = (μ + α * log(α * β)) / (1 - α)
    c3 = 1 / (1 - β)
    c4 = 1 / (1 - α * β)
    return c1 + c2 * (c3 - c4) + c4 * log(y)
end
c_star(y;p) = (1 - p.α * p.β) * y

Random.seed!(42) # for reproducible results
u(c;p) = log(c) # utility
f(k;p) = k^p.α # deterministic part of production function
OptimalGrowthModel = @with_kw (α = 0.4, β = 0.96, μ = 0.0, s = 0.1,
                  u = u, f = f, # defaults defined above
                  y = range(1e-5, 4.0, length = 200), # grid on y
                  ξ = exp.(μ .+ s * randn(250)) # monte carlo shocks
) # named tuples defaults

# test 

# True value and policy function
function v_star(y;p)
    @unpack α, μ, β = p
    c1 = log(1 - α * β) / (1 - β)
    c2 = (μ + α * log(α * β)) / (1 - α)
    c3 = 1 / (1 - β)
    c4 = 1 / (1 - α * β)
    return c1 + c2 * (c3 - c4) + c4 * log(y)
end
c_star(y;p) = (1 - p.α * p.β) * y


# test 

p = OptimalGrowthModel() # use all default parameters from named tuple
w_star = v_star.(p.y; p)  # evaluate closed form value along grid

w = T(w_star; p).w # evaluate operator, access Tw results

plt = plot(ylim = (-35,-24))
plot!(plt, p.y, w, linewidth = 2, alpha = 0.6, label = L"T(v^*)")
plot!(plt, p.y, w_star, linewidth = 2, alpha=0.6, label = L"v^*")
plot!(plt, legend = :bottomright)


# 

w = 5 * log.(p.y)  # An initial condition -- fairly arbitrary
n = 35

plot(xlim = (extrema(p.y)), ylim = (-50, 10))
lb = "initial condition"
plt = plot(p.y, w, color = :black, linewidth = 2, alpha = 0.8, label = lb)
for i in 1:n
    w = T(w; p).w
    plot!(p.y, w, color = RGBA(i/n, 0, 1 - i/n, 0.8), linewidth = 2, alpha = 0.6,
          label = "")
end

lb = "true value function"
plot!(plt, p.y, v_star.(p.y; p), color = :black, linewidth = 2, alpha = 0.8, label = lb)
plot!(plt, legend = :bottomright)

# 

function solve_optgrowth(initial_w; p, iterations = 500, m = 3, show_trace = false) 
    results = fixedpoint(w -> T(w;p).w, initial_w; iterations, m, show_trace) # Anderson iteration
    v_star = results.zero
    σ = T(results.zero;p).σ
    return (;v_star, σ, results)
end

#

initial_w = 5 * log.(p.y)
sol = solve_optgrowth(initial_w;p)
v_star_approx = sol.v_star
println("Converged in $(sol.results.iterations) to an ||residuals||_∞ of $(sol.results.residual_norm)")

plt = plot(ylim = (-35, -24))
plot!(plt, p.y, v_star_approx, linewidth = 2, alpha = 0.6,
      label = "approximate value function")
plot!(plt, p.y, v_star.(p.y;p), linewidth = 2, alpha = 0.6, label = "true value function")
plot!(plt, legend = :bottomright)

#

plt = plot(p.y, T(v_star_approx; p).σ, lw=2, alpha=0.6, label = "approximate policy function")
plot!(plt, p.y, c_star.(p.y; p), lw = 2, alpha = 0.6, label = "true policy function")
plot!(plt, legend = :bottomright)

#


############ Time iteration
using LinearAlgebra, Statistics
using BenchmarkTools, Interpolations, LaTeXStrings, Parameters, Plots, QuantEcon, Roots
using Optim, Random
using BenchmarkTools, Interpolations, Parameters, Plots, QuantEcon, Roots

function K!(Kg, g, grid, β, ∂u∂c, f, f′, shocks)
    # This function requires the container of the output value as argument Kg
    
        # Construct linear interpolation object
        g_func = LinearInterpolation(grid, g, extrapolation_bc=Line())
    
        # solve for updated consumption value
        for (i, y) in enumerate(grid)
            function h(c)
                vals = ∂u∂c.(g_func.(f(y - c) * shocks)) .* f′(y - c) .* shocks
                return ∂u∂c(c) - β * mean(vals)
            end
            Kg[i] = find_zero(h, (1e-10, y - 1e-10))
        end
        return Kg
    end
    
    # The following function does NOT require the container of the output value as argument
    K(g, grid, β, ∂u∂c, f, f′, shocks) =
        K!(similar(g), g, grid, β, ∂u∂c, f, f′, shocks)

isoelastic(c, γ) = isone(γ) ? log(c) : (c^(1 - γ) - 1) / (1 - γ)
Model = @with_kw (        α = 0.36,                            # Productivity parameter
                          β = 0.99,                            # Discount factor
                          γ = 1,                             # Risk aversion
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
Random.seed!(42) # for reproducible results.

shock_size = 250 # number of shock draws in Monte Carlo integral
shocks = collect(exp.(m.μ .+ m.s * randn(shock_size))); # generate shocks

function verify_true_policy(m, shocks, c_star)
    # compute (Kc_star)
    @unpack grid, β, ∂u∂c, f, f′ = m
    c_star_new = K(c_star, grid, β, ∂u∂c, f, f′, shocks)

    # plot c_star and Kc_star
    plot(grid, c_star, label = L"optimal policy $c^*$")
    plot!(grid, c_star_new, label = L"Kc^*")
    plot!(legend = :topleft)
end

c_star = (1 - m.α * m.β) * m.grid # true policy (c_star)
verify_true_policy(m, shocks, c_star)

function check_convergence(m, shocks, c_star, g_init; n_iter = 150)
    @unpack grid, β, ∂u∂c, f, f′ = m
    g = g_init;
    plot(m.grid, g, lw = 2, alpha = 0.36, label = L"intial condition $c(y) = y$")
    for i in 1:n_iter
        new_g = K(g, grid, β, ∂u∂c, f, f′, shocks)
        g = new_g
        plot!(grid, g, lw = 2, alpha = 0.36, label = "")
    end
    plot!(grid, c_star, color = :black, lw = 2, alpha = 0.8,
          label = L"true policy function $c^*$")
    plot!(legend = :topleft)
end

check_convergence(m, shocks, c_star, m.grid, n_iter = 15)




##

function K!(Kg, g, grid, β, ∂u∂c, f, fk′, fn', shocks)
    # This function requires the container of the output value as argument Kg
    
        # Construct linear interpolation object
        g_func = LinearInterpolation(grid, g, extrapolation_bc=Line())
    
        # solve for updated consumption value
        for (i, y) in enumerate(grid)
            function h(c)
                vals = ∂u∂c.(g_func.(f(y - c) * shocks)) .* f′(y - c) .* shocks
                return ∂u∂c(c) - β * mean(vals)
            end
            Kg[i] = find_zero(h, (1e-10, y - 1e-10))
        end
        return Kg
    end
    
    # The following function does NOT require the container of the output value as argument
    K(g, grid, β, ∂u∂c, f, fk′, fn', shocks) =
        K!(similar(g), g, grid, β, ∂u∂c, f, f′, fk′, fn', shocks)

##

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
                  f = (k,n) -> k^α * n^(1-α),          # production function
                  fk′ = k -> α * k^(α - 1),             # f′ wrt to k
                  fn′ = n -> (1-α) * (k/n)^α,             # f′ wrt to n
                  )

##

m = Model()