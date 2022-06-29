## Implementing EGM generalised to all types of shocks + finding why decision rules are flat

using Pkg
using Interpolations
using Plots
using Dolo
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
dp.grid
dp.integration_nodes
dp.integration_weights
Dolo.iweight(dp,1,1) # transition probabilities
Dolo.n_inodes(dp,1) # nb of possible future states given current state i
size_states = Dolo.n_nodes(dp) # size of matrix states


# calibration 
m_, s_, mr_, a_, x_, p_ = model.calibration[:exogenous, :states, :expectations, :poststates, :controls, :parameters]
m,s,mr,a,x,p = [SVector(e...) for e in model.calibration[:exogenous, :states, :expectations, :poststates, :controls, :parameters]]

# grid 
model.domain
Dolo.get_domain(grid)
domain = Dolo.get_domain(model)
grid.states
grid.min
grid.max 
grid = Dolo.get_discretization_options(model)
grid[:exo][:n]
grid[:endo][:n]

discretized = Dolo.discretize(model)

exogrid = discretize[1].exo
exogrid.n

endogrid = discretize[1].endo
endogrid.n

function consumption_a(model;
    dr0=Dolo.ConstantDecisionRule(model),
    φ1)
    if 
    for i in 1:
end

function egm()
    
end

