using Pkg
using Interpolations
using Plots
using Dolo: UNormal, discretize  # dependence on Dolo is actually quite minimal, we just use the exogenous shock
using StaticArrays
Pkg.add("StaticArrays")


#=
structure m qui defini le model a resoudre, methode de resolution qui ne suppose rien sur le model
# 1 var etat, 1 Control 
# let : en dors de let, 
=#

# We define the model here
m = let
    # parameter calibration

    γ = 4.0   
    σ_y = 0.1
    β = 0.96
    r = 1.038
    p = (;β, γ, r, σ_y)

    # The following computes a set of nodes / weights
    exogenous = UNormal(;σ=σ_y)
    dp = discretize(exogenous)
    x = SVector(dp.integration_nodes[:]...) # nodes 
    w = SVector(dp.integration_weights[:]...) # weights

    # Initial decision rule (consumption as a function of available wealth)
    φ0 = w -> min(w, 1.0 + 0.01*(w-1.0))

#fonction qui a w associe 

    # This is the discretization of...
    N = 10
    w_grid = range(0.8, 20; length=N) # the state-space
    a_grid = range(0.0, 20; length=N) # the post-states
    
    (;p, φ=φ0, a_grid, w_grid, integration=(w, x))

end

#=
NamedTuple 
t = (3, "JI")
length(t) # on ne peut pas modifer elements de T
 # variante NamedTuple, avec des mots cle 
 t = (a = 4, b = 3)
 t.a 
 # add LabelledArrays, se comporte comme un vector 
 # synthaxe pour creer NamedTuple 

 m.integration
 a = w-c 
 si on connait c sur la grille de ai on connait c sur la grille de wi
=#


function consumption_a(φ1, m)
    """Computes consumption at each point of the post-state grid given:
    - φ1: decision rule for consumption tomorrow
    - m:  model
    """
    
    c_a = zeros(length(m.a_grid))
    
    (weights, nodes) = m.integration
    for (n,a) in enumerate(m.a_grid)
        e = 0.0
        for i in 1:length(weights) # # βREyu'(c'(M'))
            w = weights[i] # probabilities of each nodes
            ε = nodes[i] # nodes of the stochastic process
            W = exp(ε) + a*m.p.r # M' = AR + y
            C = φ1(W) # c'(M') using c_(i-1)(.)
            e += m.p.β*w*C^(-m.p.γ)*(m.p.r) 
        end
        c_a[n] = (e)^(-1.0/m.p.γ)
    end



    return c_a

end


function egm(m; φ=m.φ, T=500, trace=false, resample=false, 
    
    τ_η=1e-8)
    """Computes the time iteration solution given:
    - φ: initial decision rule for consumption
    - m: model
    """

    logs = [] # to keep all successive decision rules

    local w_grid, c_a, φ0
    a_grid = m.a_grid

    φ0 = φ
    
    for t in 1:T
        trace ? push!(logs, deepcopy(φ0)) : nothing
        c_a = consumption_a(φ0, m) # c_new = (u')^(-1)(A)
        w_grid = a_grid + c_a # M_new = A + c_new
        c_a = min(w_grid, c_a) # c_new cannot exceed M
        φ0 = LinearInterpolation(w_grid, c_a; extrapolation_bc=Line()) # reconstructing policy function c_new(M)
    end

    if resample
        # we resample the solution so that it interpolates exactly
        # on the grid specified in m.w_grid (not the endogenous one.)
        nx = φ0(m.w_grid)
        φ0 = LinearInterpolation(m.w_grid, nx; extrapolation_bc=Line())
    end

    if trace
        return φ0, logs, w_grid, c_a
    else
        return φ0
    end

end

@time φ = egm(m; resample=true)
xvec = range(0,10;length=100)
plt = plot(xvec, xvec; xlims=(0,6))
plot!(plt, xvec, min.(xvec,φ.(xvec)))

φ = egm(m; resample=false)
φs = egm(m; resample=true)
plot(xvec, xvec; label="w")
plot!(φ.itp.knots[1], φ.itp.coefs; marker= "o", label="c(W) endogenous")
plot!(φs.itp.ranges[1], φs.itp.itp.coefs; marker= "o", label="c(A) fixed", xlabel="State w", ylabel="Control c(w)")


# plotting 
@time soltrace = egm(m; resample=true, trace = true) 
@time sol = egm(m; resample=true)
xvec = range(0,10;length=100)
function convergenceEGM(soltrace, sol)
        soltrace = soltrace
        trace = soltrace[2]
        sol = sol
        plt = plot()
        plot!(plt, xvec, xvec; label="w", ylims=(0,1.5))
    for i=1:20:length(trace)
        plot!(plt, sol.itp.ranges[1], trace[i]; marker= "o")
    end
    plot!(plt, xlabel = "Wealth", ylabel = "Consumption")
    plot!(plt, legend = false)
end
convergenceEGM(soltrace, sol)

# replace iid shock with markov chain

