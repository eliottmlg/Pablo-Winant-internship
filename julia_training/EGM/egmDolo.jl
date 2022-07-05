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


# declaring model
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/Dolo.jl/examples/models/consumption_savings_iid.yaml"
readlines(filename)
model = yaml_import(filename)
model.symbols

#equations
F_tran = Dolo.get_factory(model, "transition")
F_arb = Dolo.get_factory(model, "arbitrage")[1]
F_g = Dolo.get_factory(model, "half_transition")
F_τ = Dolo.get_factory(model, "direct_response_egm")
F_h = Dolo.get_factory(model, "expectation")
F_aτ =  Dolo.get_factory(model, "reverse_state")

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
# N_m = max(Dolo.n_nodes(grid_exo),1) will not work for IID processes because n_nodes(grid_exo) not defined
if typeof(dp) != Dolo.DiscretizedIIDProcess
    size_states = Dolo.n_nodes(dp)
else 
    size_states = size(dp.integration_nodes,1)
end

# initial policy function
#φ0 = Vector{Any}(undef, size_states)
function φfunction(i, ss::SVector{1})::SVector # i is the current exogenous state, takes state at i
    xx = min(ss[1], 1.0 + 0.01*(ss[1]-1.0))
    return ss, xx
end

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

function consumption_a(model,φ)
    φ1 = φ
    c_a = Matrix{typeof(x)}(undef, length(a0), size_states)
    
    for i in 1:size_states
        for (n,a) in enumerate(a0) 
            
            m = SVector{length(Dolo.node(dp,i))}(Dolo.node(dp,i)) # convert to Svector{Float64}
            rhs = mr*0.0
            zz = zeros(size_states)
                        
            for j in 1:size_states
                
                #=
                if typeof(dp) != Dolo.DiscretizedIIDProcess
                    m = Dolo.inode(dp)[j,1]
                    w = Dolo.iweight(dp,)
                else 
                    x = dprocess.integration_nodes[j,1]
                    w = dprocess.integration_weights[j]

                end
                =#

                #inc = states[j,1]
                #prob = transitions[i,j]
                M = SVector{length(Dolo.inode(dp,i,j))}(Dolo.inode(dp,i,j)) # convert to Svector{Float64}
                w = Dolo.iweight(dp,i,j)

                #W = exp(inc) + a*m.p.r # M' = AR + y
                #ss = g(m,a,M,p) # S = g(m,a,M), half_transition   NEEDS 5 arguments
                ss = g(m,a,M,p) # S = g(m,a,M), half_transition

                # xx = φ1[j](ss) # c'(M') using c_(i-1)(.)  
                ss, xx = φ1(j,ss) # c'(M') using c_(i-1)(.)  
                
                #Φ[j] = m.p.β * prob * (C/m.p.cbar)^(-m.p.γ) * (m.p.r) 
                zz[j] = w * h(M,ss,xx,p) # z = E(h(M,S,X)), expectation

            end
            #rhs = LinearAlgebra.dot(Φ, transitions[i,:])
            rhs += sum(Dolo.iweight(dp,i,j) * zz[j] for j=1:size_states)

            #c_a[n, i] = m.p.cbar * (rhs)^(-1.0/m.p.γ)
            τ(m,a,zz,p) #print
            c_a[n,i] = τ(m,a,zz,p) # c = \tau(m,a,z), direct_response_egm
        end    
    end
end

consumption_a(model, φfunction)

# interpolation 
itp = Vector{Any}(undef,size_states) # empty vector to host interpolation objects 

struct MyDR # creating a structure for neater incorporation
    itp::Vector{Any}
end
(mydr::MyDR)(i,s) = mydr.itp[i](s)

function egm(model; φfunction=nothing, T=500, trace=false, resample=false, τ_η=1e-8)
    #states = m.markovpricess[1]
    inodes_vec = dprocess.integration_nodes
    logs = []
    #local w_grid, c_a, φ0
    #w_grid = s0
    local c_a, φ0

    φ0 = φfunction

    for t in 1:T
        trace ? push!(logs,deepcopy(itp)) : nothing
        c_a = consumption_a(model,φ0) # c_new = (u')^(-1)(A)

        for i in 1:length(inodes_vec) 
            #w_grid = a_grid + c_a[:,i] # M_new = A + c_new
            w_grid = aτ(m,a,c_a[:,i],p) # s = a\tau(m,a,x), reverse_state
            c_a[:,i] = min.(w_grid, c_a[:,i]) # c_new cannot exceed M
            itp[i] = LinearInterpolation(w_grid, c_a[:,1], extrapolation_bc = Line())
            #φ0[i] = LinearInterpolation(w_grid, c_a[:,i]; extrapolation_bc=Line()) # reconstructing policy function c_new(M)
        end
        res = MyDR(itp)
        mydr(i,w_grid) = res.itp[i](w_grid)
        # mydr(1,x)
    end
        if resample
            for i in 1:size_states
            # we resample the solution so that it interpolates exactly
            # on the grid specified in m.w_grid (not the endogenous one.)
            nx = min.(w_grid, mydr(i,w_grid))
            itp[i] = LinearInterpolation(w_grid, nx, extrapolation_bc=Line())
            end
        end
        res = MyDR(itp)
        mydr(i,w_grid) = res.itp[i](w_grid)
        # mydr(1,x)

        if trace
            return res, mydr, logs, w_grid, c_a
        else
            return res, mydr # 3-element vector, each element is an object of linear interpolation containing grid and values on the grid (10-element vector containing optimal dr(i,a))
        end
end

@time φs = egm(model, φ)


