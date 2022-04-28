
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
using LinearAlgebra, Statistics
using Plots, QuantEcon, Interpolations, NLsolve, Optim, Random, Parameters
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