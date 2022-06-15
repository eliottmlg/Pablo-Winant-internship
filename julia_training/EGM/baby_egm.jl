## Implementing EGM generalised to all types of shocks + finding why decision rules are flat

using Pkg
using Interpolations
using Plots
using Dolo
using StaticArrays
using QuantEcon
using LinearAlgebra

# We define the model here
m = let
    # parameter calibration

    γ = 4.0   
    σ_y = 0.01
    β = 0.96
    r = 1.038
    rho = 0.95
    N_mc = 3
    sigma = Array{Float64}(undef, 1, 1)
    sigma[1,1] = σ_y
    p = (;β, γ, r, σ_y, rho)

    # Markov process
    ar1 = Dolo.VAR1(ρ = rho, Σ = sigma)
    mc = Dolo.discretize(ar1; n=N_mc)
    states = mc.values
    transitions = mc.transitions

    # Initial decision rule (consumption as a function of available wealth)
    φ0 = Vector{Any}(undef, N_mc)
    for i in 1:size(states,1)
    φ0[i] = w -> min(w, 1.0 + 0.01*(w-1.0))
    end
    φ0  # vector of length 3 of anonymous functions

    # This is the discretization of...
    N = 100
    w_grid = range(0.8, 20; length=N) # the state-space
    a_grid = range(0.0, 20; length=N) # the post-states
    
    (;p, φ=φ0, a_grid, w_grid, markovprocess=(states, transitions))

end


function consumption_a(φ1, m)
    """Computes consumption at each point of the post-state grid given:
    - φ1: decision rule for consumption tomorrow
    - m:  model
    """

    (states, transitions) = m.markovprocess
    c_a = Matrix{Float64}(undef, length(m.a_grid), size(states,1))

    for i in 1:size(states,1)  # βREyu'(c'(M'))

        for (n,a) in enumerate(m.a_grid)

            rhs = zeros(size(states,1))
            Φ = zeros(size(states,1))

            for j in 1:size(states,1)

                inc = states[j,1]
                prob = transitions[i,j]
                W = exp(inc) + a*m.p.r # M' = AR + y
                C = φ1[j](W) # c'(M') using c_(i-1)(.)  
                Φ[j] = m.p.β * prob * C^(-m.p.γ)*(m.p.r) 

            end

            rhs = LinearAlgebra.dot(Φ, transitions[i,:])
            c_a[n, i] = (rhs)^(-1.0/m.p.γ)
        end 
    end
    return c_a # matrix length(a_grid) x length(states) containing updated consumption levels (Float64)
end


function egm(m; φ=m.φ, T=500, trace=false, resample=false, τ_η=1e-8)
    """Computes the time iteration solution given:
    - φ: initial decision rule for consumption
    - m: model
    """

    states = m.markovprocess[1]
    logs = [] # to keep all successive decision rules

    local w_grid, c_a, φ0

    a_grid = m.a_grid
     φ0 = φ
    
    for t in 1:T
        
        trace ? push!(logs, deepcopy(φ0)) : nothing
        c_a = consumption_a(φ0, m) # c_new = (u')^(-1)(A)
        
        for i in 1:size(states,1)
                
            w_grid = a_grid + c_a[:,i] # M_new = A + c_new
            c_a[:,i] = min(w_grid, c_a[:,i]) # c_new cannot exceed M
            φ0[i] = LinearInterpolation(w_grid, c_a[:,i]; extrapolation_bc=Line()) # reconstructing policy function c_new(M)
        end
    end

    if resample
        for i in 1:size(states,1)
            # we resample the solution so that it interpolates exactly
            # on the grid specified in m.w_grid (not the endogenous one.)
            nx = φ0[i](m.w_grid)
            φ0[i] = LinearInterpolation(m.w_grid, nx; extrapolation_bc=Line())
        end
    end

    if trace
        return φ0, logs, w_grid, c_a
    else
        return φ0 # 3-element vector, each element is an object of linear interpolation containing grid and values on the grid (10-element vector containing optimal dr(i,a))
    end

end

@time φs = egm(m; resample=true)
function result(φs)
    xvec = range(0,10;length=10)
    plt = plot(xvec, xvec; xlims=(0,10), xlabel="State w", ylabel="Control c(w)")
    for i in 1:length(m.markovprocess[1])
        x = φs[i].itp.ranges[1]
        plt = plot!(φs[i].itp.ranges[1], min.(x,φs[i](x)); marker= "o")
    end
    plt
end 
result(φs)


