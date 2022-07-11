## Implementing EGM generalised to all types of shocks + finding why decision rules are flat

using Pkg
using Interpolations
using Plots
using Dolo
using Dolo: UNormal, discretize 
using StaticArrays
using QuantEcon
using LinearAlgebra
import Dolang

# Changes to Dolo.jl, so need to run the following:
# model.jl: 
# discretizemodel()
# get_factory_noloopEGM()

# declaring model
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/Dolo.jl/examples/models/consumption_savings_iid.yaml"
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/julia_training/EGM/consumption_savings_ar1.yaml"

function DeclareModel(filename)
    readlines(filename)
    model = yaml_import(filename)    
end
DeclareModel(filename)

#equations
function ModelParse(model)

    #equations
    F_tran = Dolo.get_factory_noloopEGM(model, "transition")
    F_arb = Dolo.get_factory_noloopEGM(model, "arbitrage")[1]
    F_g = Dolo.get_factory_noloopEGM(model, "half_transition")
    F_τ = Dolo.get_factory_noloopEGM(model, "direct_response_egm")
    F_h = Dolo.get_factory_noloopEGM(model, "expectation")
    F_aτ =  Dolo.get_factory_noloopEGM(model, "reverse_state")

    code_tran = Dolang.gen_generated_gufun(F_tran)
    tran = eval(code_tran)
    code_arb = Dolang.gen_generated_gufun(F_arb)
    arbi = eval(code_arb)
    code_g = Dolang.gen_generated_gufun(F_g)
    g = eval(code_g)
    code_τ = Dolang.gen_generated_gufun(F_τ)
    τ = eval(code_τ)
    code_aτ = Dolang.gen_generated_gufun(F_aτ)
    aτ = eval(code_aτ)
    code_h = Dolang.gen_generated_gufun(F_h)
    h = eval(code_h)

    
    # exogenous shock
    shock = Dolo.get_exogenous(model)
    dp = Dolo.discretize(shock)
    size_states = Dolo.n_inodes(dp,1)

        
    # calibration 
    m_, s_, mr_, a_, x_, p_ = model.calibration[:exogenous, :states, :expectations, :poststates, :controls, :parameters]
    m,s,mr,a,x,p = [SVector(e...) for e in model.calibration[:exogenous, :states, :expectations, :poststates, :controls, :parameters]]

    # grids
    grid, dprocess = Dolo.discretize(model)
    grid_endo = grid.endo
    grid_exo = grid.exo
    grid_fixed = grid_endo
    s0 = Dolo.nodes(grid_endo)
    a0 = Dolo.nodes(grid_fixed)

end
ModelParse(model)

Dolo.reverse_state(model, m,a,m,p)

# initial policy function
function φfunction(i, ss::SVector{1}, itp)::SVector # i is the current exogenous state, takes state at i
    #xx = min(ss[1], 1.0 + 0.01*(ss[1]-1.0))
    xx = 0.9*ss[1]
    xx = SVector(xx...)
    return ss, xx
end # to call: φfunction(1,s0[1])

function consumption_a(model,φ=nothing)
    φ1 = φ
    c_a = Matrix{typeof(x)}(undef, length(a0), size_states)
    
    for i in 1:size_states
        for (n,a) in enumerate(a0) 
            
            m = SVector(Dolo.inode(dp,i,i)...)
            zz = mr*0.0
            zzj = zeros(SVector{1,Float64}, size_states)
                        
            for j in 1:size_states
                
                M = SVector(Dolo.inode(dp,i,j)...)
                w = Dolo.iweight(dp,i,j)

                ss = g(m,a,M,p) # S = g(m,a,M), half_transition

                ss, xx = φ1(i,ss) # c'(M') using c_(i-1)(.)  
                
                zzj[j] = w * h(M,ss,xx,p) # z = E(h(M,S,X)), expectation

            end
            zz += sum(zzj[j] for j=1:size_states)
            τ(m,a,zz,p)
            c_a[n,i] = τ(m,a,zz,p) # c = \tau(m,a,z), direct_response_egm
        end    
    end
    return c_a
end

# computing updated deicisions after from initial decision rule to iteration one  
cprime = consumption_a(model, φfunction)# same consumption for all state of m
# expected since shock is iid

# iteration by iteration

# t = 1
cprime = consumption_a(model,φfunction) # c_new = (u')^(-1)(A)
w_grid = s0
        for i in 1:size_states
            for n in 1:length(w_grid)
                w_grid[n] = aτ(m,a,cprime[n,i],p) # s = a\tau(m,a,x), reverse_state
                cprime[n,i] = min(w_grid[n], cprime[n,i]) # c_new cannot exceed M
            end
            cprimeFloat = reinterpret(Float64, cprime)
            s0Float = reinterpret(Float64,w_grid)
            itp[i] = LinearInterpolation(s0Float, cprimeFloat[:,i], extrapolation_bc = Line())
        end
cprime # exactly the same because
itp[1]
s0Float
function plotcprime1(itp)
    plt = plot(ylims = 1.2)
    plot!(plt, s0Float,s0Float; legend = true)
    plot!(plt, s0Float, itp[1]; marker = "o")
    plt
end
plotcprime1(itp)

# t = 2

cprime = consumption_a(model,itp) # c_new = (u')^(-1)(A)
w_grid = s0
        for i in 1:size_states
            for n in 1:length(w_grid)
                w_grid[n] = aτ(m,a,cprime[n,i],p) # s = a\tau(m,a,x), reverse_state
                cprime[n,i] = min(w_grid[n], cprime[n,i]) # c_new cannot exceed M
            end
            cprimeFloat = reinterpret(Float64, cprime)
            s0Float = reinterpret(Float64,w_grid)
            itp[i] = LinearInterpolation(s0Float, cprimeFloat[:,i], extrapolation_bc = Line())
        end
cprime # exactly the same because
itp[1]
s0Float
function plotcprime1(itp)
    plt = plot(ylims = 1.2)
    plot!(plt, s0Float,s0Float; legend = true)
    plot!(plt, s0Float, itp[1]; marker = "o")
    plt
end
plotcprime1(itp)



# plotting for exogenous state 1 updated consumption along fixed grid A 
function plotcprime1(cprime)
    plt = plot()
    for i in 1:size_states
        plot!(plt, a0, cprime[:,i]; marker = "o", legend = false)
    end
    plt
end
plotcprime1(cprime)

# interpolation 
itp = Vector{Any}(undef,size_states) # empty vector to host interpolation objects 
s0Float = reinterpret(Float64, s0) # convert from Vector(Svector(Float64)) to Vector(Float64)
cprimeFloat = reinterpret(Float64, cprime) # for function LinearInterpolation

struct MyDR # creating a structure for neater incorporation
    itp::Vector{Any}
end
(mydr::MyDR)(i,s) = mydr.itp[i](s)

# interpolated policy functions
mydr(i,s) = φs[2].itp[i](s)
mydrlogs(i,s,t) = φs[3][t].itp[i](s)

function toyegm(model; φfunction=nothing, T=500, trace=false, resample=false, τ_η=1e-8)
    logs = []
    
    local cprime, φ0

    φ0 = φfunction
    w_grid = s0

    for t in 1:T
        cprime = consumption_a(model,φ0) # c_new = (u')^(-1)(A)
        for i in 1:size_states
            for n in 1:length(w_grid)
                w_grid[n] = aτ(m,a,cprime[n,i],p) # s = a\tau(m,a,x), reverse_state
                cprime[n,i] = min(w_grid[n], cprime[n,i]) # c_new cannot exceed M
            end
            cprimeFloat = reinterpret(Float64, cprime)
            s0Float = reinterpret(Float64,w_grid)
            itp[i] = LinearInterpolation(s0Float, cprimeFloat[:,i], extrapolation_bc = Line())
        end
        trace ? res = MyDR(itp) : nothing
        trace ? push!(logs,deepcopy(res)) : nothing
    end
    dr = cprime
    res = MyDR(itp)
    
    if trace
        return dr, res, logs
    else
        return dr, res
    end     
end

# policy functions for each state
@time φs = toyegm(model; φfunction)
mydr(1,s0Float)
function policyfunction(size_states)
    constraint = (1:length(s0Float))
    plt = plot()
    #plot!(plt, constraint, constraint, xlims=s0Float[1000])
    for i in 1:size_states
        plot
        plot!(plt, s0Float, mydr(i,s0Float); legend = true)
    end
    plt
end
policyfunction(size_states)


# policy functions at each iteration
@time φs = toyegm(model; φfunction, trace = true)
mydrlogs(1,s0Float,500)
function iterationEGM(T) 
    plt = plot()
    for t in 1:20:T
        plot!(plt, s0Float, mydrlogs(1,s0Float,t); legend = true)
    end
    plt
end
iterationEGM(500)

