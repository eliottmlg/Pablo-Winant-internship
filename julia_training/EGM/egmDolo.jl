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

Dolo.get_factory(model, "half_transition")
specs = RECIPES[:dtcc][:specs]
Symbol(half_transition)
keys(specs)

filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/Dolo.jl/examples/models/rbc_mc.yaml"
readlines(filename)
model = yaml_import(filename)
model.symbols
specs = RECIPES[:dtcc][:specs]
specs[:half_transition]
recipe = specs[Symbol("half_transition")]
this = recipe[:eqs]
enumerate(this)
this[1][3]

arguments = OrderedDict(
    Symbol(l[3]) => [stringify(e,l[2]) for e in symbols[Symbol(l[1])]]
    for l in recipe[:eqs] if !(l[1]=="parameters")
)


collect(keys(F_g.arguments))
symbols = Dolo.get_symbols(model)
symbols[Symbol(:poststates)]


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
N_m = max(Dolo.n_nodes(grid_exo),1)
if typeof(dp) != Dolo.DiscretizedIIDProcess
    size_states = Dolo.n_nodes(dp)
else 
    size_states = size(dp.integration_nodes,1)
end

# initial policy function
φ0 = Vector{Any}(undef, size_states)
function φ1(i, s::SVector)::SVector # i is the current exogenous state, takes state at i
    #s::[w] returns these irrespective of i
    #x::[c] 
    return s, x
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

# recursively make all keys at any layer of nesting a symbol
# included here instead of util.jl so we can call it on RECIPES below
_symbol_dict(x) = x
_symbol_dict(d::AbstractDict) =
    Dict{Symbol,Any}([(Symbol(k), _symbol_dict(v)) for (k, v) in d])
const src_path = dirname(@__FILE__)
const pkg_path = dirname(src_path)
Pkg.add("YAML")
import YAML; using YAML: load_file, load
const RECIPES = _symbol_dict(load_file(joinpath(src_path, "recipes.yaml")))
specs = RECIPES[:dtcc][:specs]

function consumption_a(model,φ1)
    φ1 = φ1
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

consumption_a(model, φ1)


function egm(model; φ=nothing, T=500, trace=false, resample=false, τ_η=1e-8)
    #states = m.markovpricess[1]
    inodes_vec = dprocess.integration_nodes
    logs = []
    local w_grid, c_a, φ0
    a_grid = a0
    φ0 = φ

    for t in 1:T
        trace ? push!(logs,deepcopy(φ0)) : nothing
        c_a = consumption_a(model,φ0) # c_new = (u')^(-1)(A)

        for i in 1:length(inodes_vec) 
            #w_grid = a_grid + c_a[:,i] # M_new = A + c_new
            w_grid = aτ(m,a,c_a[:,i],p) # s = a\tau(m,a,x), reverse_state
            c_a[:,i] = min(w_grid, c_a[:,i]) # c_new cannot exceed M
            φ0[i] = LinearInterpolation(w_grid, c_a[:,i]; extrapolation_bc=Line()) # reconstructing policy function c_new(M)
        end
    end
end

struct MyDR
    itp::Vector{Any}
end

(mydr::MyDR)(i,s) = mydr.itp[i](s)
SVector(1...)
SVector(Dolo.node( dp,2)...)
mr
zeros(typeof(mr), 10)

function φ(i, s::SVector)::SVector # i is the current exogenous state, takes state at i
    #s::[w] returns these irrespective of i
    #x::[c] 
    return s, x
end
φ(1,s0[1])

SVector(1,2,3)
SVector{2,Any}(x1,x2,...)
consumption_a(φ)   #Vector{SVector}#

function φ(i, s::Vector{T})::SVector{T} where T
    # s: [w]
    # x: [c]
    
end

@time φs = egm(model, φ)

