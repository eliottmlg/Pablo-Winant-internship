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
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/julia_training/consumption-savings/CS_py.yaml"
readlines(filename)
model = yaml_import(filename)

#equations
arb = Dolo.get_equation_block(model, "arbitrage")
tran = Dolo.get_assignment_block(model, "transition")
h = Dolo.get_assignment_block(model, "expectation")
g = Dolo.get_assignment_block(model, "half_transition")
τ = Dolo.get_assignment_block(model, "direct_response_egm")
aτ = Dolo.get_assignment_block(model, "reverse_state")


F_tran = Dolo.get_factory(model, "transition")
F_arb = Dolo.get_factory(model, "arbitrage")[1]
F_g = Dolo.get_factory(model, "half_transition")
F_τ = Dolo.get_factory(model, "direct_response_egm")
F_h = Dolo.get_factory(model, "expectation")
F_aτ =  Dolo.get_factory(model, "reverse_state")

vector_F = [F_arb, F_tran, F_g, F_τ, F_aτ, F_h]
code = Dolang.gen_generated_gufun(F_tran)
tran = eval(code)
code = Dolang.gen_generated_gufun(F_arb)
arbi = eval(code)
code = Dolang.gen_generated_gufun(F_g)
g = eval(code)
code = Dolang.gen_generated_gufun(F_τ)
τ = eval(code)
code = Dolang.gen_generated_gufun(F_aτ)
aτ = eval(code)
code = Dolang.gen_generated_gufun(F_h)
h = eval(code)

# exogenous shock
shock = Dolo.get_exogenous(model)
dp = Dolo.discretize(shock, n = grid[:exo][:n])
N_m = max(Dolo.n_nodes(grid_exo),1)
φ0 = Vector{Any}(undef, N_m)



typeof(dp)
dp.grid
dp.integration_nodes
dp.integration_weights
Dolo.iweight(dp,1,1) # transition probabilities
Dolo.n_inodes(dp,1) # nb of possible future states given current state i
size_states = Dolo.n_nodes(dp) # size of matrix states


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

model.domain
Dolo.get_domain(grid)
domain = Dolo.get_domain(model)
grid.states
grid.min
grid.max 
grid = Dolo.get_discretization_options(model)
grid[:exo][:n]
grid[:endo][:n]

grid, dproces = Dolo.discretize(model)
grid_endo = grid.endo
grid_exo = grid.exo



function consumption_a(model;
    dr0=Dolo.ConstantDecisionRule(model),
    φ1)

    #nodes = Dolo.inode(dprocess,1)

    if typeof(dp) != Dolo.DiscretizedIIDProcess
        size_states = Dolo.n_nodes(dp)
    else 
        size_states = size(dp.integration_nodes,1)
    end
    c_a = Matrix{Float64}(undef, length(a_grid), size_states)

    for i in 1:size_states
        for (n,a) in enumerate(a0) 
        
            rhs = zeros(size_states)
            Φ = zeros(size_states)

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
                inode = Dolo.inode(dp,i,j) 
                weight = Dolo.iweight(dp,i,j)

                #W = exp(inc) + a*m.p.r # M' = AR + y
                state = g(inode,s,x,inode,p) # s = g(m,a,m,p), half_transition

                C = φ1[j](W) # c'(M') using c_(i-1)(.)  
                Φ[j] = m.p.β * prob * (C/m.p.cbar)^(-m.p.γ) * (m.p.r) 

            end
        end
        
    end
end

function egm()
    
end

