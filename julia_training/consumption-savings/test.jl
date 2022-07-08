
struct q; w; e; end
r = q(4,"Hi")
@unpack w, e = r
w, e

struct A; a; b; c; end
d = A(4,7.0,"Hi")
@unpack a, c = d

grid_max = 16
grid_size = 50
b = 0.0
asset_grid = range(-b, grid_max, length = grid_size)

x = (rand(10))
y = (randn(10))
longueur = 1:length(x)

polynomial(x) = x^2 * 0.5 + x * 5 + 1 
polynomial(10)
plot(1:10, polynomial.(y), label = "polynomial of y")
plot!(1:10, y, label = "y")

similar(x)

# ex1 

# Setup
cp = ConsumerProblem()
N = 100

# VI
V, c = initialize(cp)
println("Starting value function iteration")
V
c
c[25,1]
plot(cp.asset_grid, c[:,1])
consumptionVI = (rand(N)) 

for i in 1:N
    V = T(cp, V)
    consumptionVI[i] = T(cp, V, ret_policy=true)[25,1]
end

## why VI consumption drops from 8.41 to 0.75 suddenly 
V, c = initialize(cp)
c[25,1]
V = T(cp, V)
c_afteroneiteration = T(cp, V, ret_policy=true)[25,1]
##

c1 = T(cp, V, ret_policy=true)
consumptionVI

# TI
V2, c2 = initialize(cp)
println("Starting policy function iteration")
c2
c2[25,1]
consumption45line = c2
plot(cp.asset_grid, c2[:,1])
consumptionTI = (rand(N)) 
K(cp,c2)

for i in 1:N
    c2 = K(cp, c2)
    consumptionTI[i] = c2[25,1]
end

K(cp, c2)
consumptionTI

#plotting convergence of policy function 
plot(1:N, consumptionTI, label = "TI consumption low income")
plot!(1:N, consumptionVI, label = "VI consumption low income")
plot!(xlabel = "number of iterations", ylabel = "Consumption level")
plot!(legend = :topright)




##### time iteration Algorithm

function time_iteration(model;
    dr0=Dolo.ConstantDecisionRule(model),
    discretization=Dict(),
    interpolation=:cubic,
    verbose=true,
    details=true,
    ignore_constraints=false,
    trace = false,
    tol_η = 1e-8,
    tol_ε = 1e-8,
    maxit=500
)
    
    F = Euler(model; discretization=discretization, interpolation=interpolation, dr0=dr0,  ignore_constraints=ignore_constraints)

    complementarities = false


    if interpolation != :cubic
        error("Interpolation option ($interpolation) is currently not recognized.")
    end

    ti_trace = trace ? IterationTrace([x0]) : nothing


    z0 = deepcopy(F.x0)
    z1 = deepcopy(z0)

    log = IterationLog(
        it = ("n", Int),
        err =  ("ϵₙ=|F(xₙ,xₙ)|", Float64),
        sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
        lam = ("λₙ=ηₙ/ηₙ₋₁",Float64),
        elapsed = ("Time", Float64)
    )

    initialize(log, verbose=verbose; message="Time Iteration")

    local err_η, z0, z1

    err_η_0 = NaN
   # err_η = err_η_0 #
   # err_ε_0 = NaN #
   # err_ε = err_ε_0 #

    it = 0

    while it<=maxit

        it += 1

        t1 = time_ns()

        r0 = F(z0, z0; set_future=true)

        err_ε = norm(r0)
       # if err_ε<tol_ε
       #    break
       # end

        fun = u->F(u, z0; set_future=false)
        dfun = u->df_A(F,u, z0; set_future=false)

        sol = newton2(
            fun, dfun,
            z0, verbose=false
        )

        z1 = sol
        δ = z1 - z0

        trace && push!(ti_trace.trace, z1)

        err_η = norm(δ)
        gain = err_η / err_η_0
        err_η_0 = err_η



        # z0.data[:] .= z1.data
        z0 = z1

        elapsed = time_ns() - t1

        elapsed /= 1000000000

        append!(log; verbose=verbose, it=it, err=err_ε, sa=err_η, lam=gain, elapsed=elapsed)


        if err_η<tol_η && err_ε<tol_ε
            break 
        end
        
    end

    finalize(log, verbose=verbose)


    res = TimeIterationResult(F.dr.dr, it, complementarities, F.dprocess, err_η<tol_η, tol_η, err_η, log, ti_trace)
    return res

end

# iteration 110 is not what is evaluated by the convergence dummy
# the 111th iteration would be evaluated, but does not show up in the table
#  SO, either need to change order in the TI function
#  or just change the displaying


for i in 1:2
    if i==2
        x=1
        print(x)
    end
    if i==1
        x=0
        print(x)
    end
end     # variable x requires a new binding each time a loop is executed 



function volume_sphere(r)
   # (4/3)(round(π, sigdigits=6))(r)^3
    (4/3) * round(π, sigdigits=6) * r^3
  end

  function volume_sphere(r)
    (4/3)*π*r^3
  end

println(volume_sphere(5))





x = 3.141592653589793238462643383279502884197
xmax = maximum(x)-30
xmaxconverted = trunc(Int32, xmax) 
print(typeof(xmaxconverted)) 







x = 1.659259375937593759357

floor(x)

y = x - floor(x)

x = round(Int, x)
x


vector_time = rand(500)
vector_time[1] = 5
vector_time

heurenow = time_ns() 
duree = time_ns()-heurenow
duree /= 1000000000

x=0.5
if x>=1
    y = NaN
elseif x<1
    y = 4*x 
end

a=4
b=5
function compare(a, b)
    a == b && return "equal to"
    a < b ? "less than" : "greater than"
end

a == b && return "equal to"
a < b ? "less than" : "greater than"

ti_trace = trace ? IterationTrace([F.x0]) : nothing
indeed = false
a = indeed ? 3 : nothing

c = Array{Union{Float64,Nothing},2}(nothing, 3, 2)
c


mutable struct concrete
    pillar1::Array{Any,1}
end
mutable struct bridge
    matrice1::Union{Nothing,concrete}
end

vector1 = [1, 1, 1]
typeof(vector1)
example = bridge(concrete(vector1))
example
example.matrice1.pillar1 = [1, 2, 3]
example.matrice1.pillar1

concrete(vector1)

true && push!(example.matrice1.pillar1, 1) #
pop!(example.matrice1.pillar1)

vector1 = [1, 1, 1]
length(vector1)

beton = concrete(vector1)
return beton



# AxisArray

V = AxisArray(rand(10); row='a':'j')  # AxisArray(rand(10), Axis{:row}('a':'j'))
V[row='c'] == V[Axis{:row}('c')] == V[row=3] == V[3]
AxisArrays.axes([1,2])
rand(10)

using Pkg; pkg"add AxisArrays Unitful"
using AxisArrays, Unitful, Random
fs = 40000;
import Unitful: s, ms, µs
rng = Random.MersenneTwister(123)
y = randn(rng, 60*fs+1)*3

for spk = (sin.(0.8:0.2:8.6) .* [0:0.01:.1; .15:.1:.95; 1:-.05:.05] .* 50,
    sin.(0.8:0.4:8.6) .* [0:0.02:.1; .15:.1:1; 1:-.2:.1] .* 50)
i = rand(rng, round(Int,.001fs):1fs)
while i+length(spk)-1 < length(y)
 y[i:i+length(spk)-1] += spk
 i += rand(rng, round(Int,.001fs):1fs)
end
end

A = AxisArray(hcat(y, 2 .* y); time = (0s:1s/fs:60s), chan = ([:c1, :c2]))
hcat(y, 2 .* y)
(0s:1s/fs:60s)
([:c1, :c2])
petit = 1/fs
60/petit
A[time=1]
A[chan = :c2, time = 1:5]


model.symbol
show(stdout, "text/plain", model.symbols)
sol = Dolo.time_iteration(model; trace = true)
sol

matrice = [1 2 3; 4 5 6]
matrice[:,2]


# 
function plottingDR(res)
    grid = length(res.ti_trace)
    DR = res.ti_trace

    plt = plot()
    plot!(plt, grid, DR, label = "decision rule")
    plot!(plt, legend = :topleft)
end

# pull up request, loop for plotting, read paper, read endogenous grid points, 

hello = 0.1
random = UNormal(;σ=hello)
new = discretize(random)
uncertain = SVector(new.integration_nodes[:]...)

using Dolo: UNormal, discretize  # dependence on Dolo is actually quite minimal, we just use the exogenous shock
σ_y = 0.1
exogenous = UNormal(;σ=σ_y)
exogenous
typeof(UNormal(;σ=σ_y))
dp = discretize(exogenous, n=5)
x = SVector(dp.integration_nodes[:]...) # nodes 
w = SVector(dp.integration_weights[:]...) # weights
integration=(w, x)
typeof(integration)
length(integration)

φ = w -> min(w, 1.0 + 0.01*(w-1.0))
φ1 = φ  # it is a anonymous function, not a vector
C = φ1(3)

φ0 = Vector{Any}(undef, 3)
φ0
for i in 1:3
φ0[i] = w -> min(w, 1.0 + 0.01*(w-1.0))
end 
φ0
φ0 = reduce(hcat, φ0)
φ0[1,1](4)

lol = w -> min(w, 1.0 + 0.01*(w-1.0))

a_grid = range(0.0, 20; length=N=10)
c_a = zeros(length(m.a_grid))
(weights, nodes) = m.integration
for (n,a) in enumerate(m.a_grid)
    return n, a 
end


# stepping through pairs of from two sequences
countries = ("Japan", "Korea", "China")
cities = ("Tokyo", "Seoul", "Beijing")

for (country, city) in zip(countries, cities)
    println("The capital of $country is $city")
end

for (i, country) in enumerate(countries) # (index 1, country 1)
    city = cities[i]
    println("The capital of $country is $city")
end

vide = []

x = 1
either = x -> min(x, 1.0 + 0.01*(x-1.0))
either # anonymous function named either 





a = [[1,2], [2,3], [3,4]]
reduce(hcat, a)
z = zeros(2)
zeroo = [z, z]
reduce(hcat, zeroo)
zeros(10)
Matrix{Float64}(undef, 10, 10)

hello = Vector{Vector{Float64}}(undef, 3)
for i in 1:3
    hello[i] = [1,2,3]
end
hello
ca = reduce(hcat, hello)


c_a = consumption_a(φ0, m) # c_new = (u')^(-1)(A)
for i in 1:3 
    w_grid = a_grid + c_a[:,i] # M_new = A + c_new
    c_a[:,i] = min(w_grid, c_a[:,i]) # c_new cannot exceed M
    φ0 = LinearInterpolation(w_grid, c_a; extrapolation_bc=Line())   
end


@time φt = egm(m; resample=true, trace = true)
φt[1].itp.knots # error this obejct does not exist with trace 
φt[1][1].itp.ranges # ranges with resample
@time φs = egm(m; resample=true)
φs[1].itp.ranges
φs[2].itp.ranges[1]
φs[9].itp.ranges[1]
(20-0.8)/.1939 # approx 100, length of grid
φs
@time φ = egm(m; resample=false)
φ[1].itp.knots[1]


φt[2][1][5][1]
length(φt[2])
φt[2]
φt[2][1][1]
φt[2][500][1]


length(φt)
φt.w_grid
φt[1] # 3-element vector, for each markov states, each element is a 50-element vector, consumption decision on the grid of 50 points
φt[2] # 500-element vector, for 500 iterations, each element is a 50-element vector
φt[3] # 50-element vector, each element is a number, w_grid at the 500th iteration!
φt[2][500] # 3-element vector, last iteration, for each markov states, consumption along the grid 
φt[2][500][1]
φt[4] # 50x3 matrix, c_a at the 500th iteration
φt[4][:,1]
φt[1][1]-φt[4][:,1] # why not the same ?


x = φs[1].itp.ranges[1]
xvec = range(0,m.w_grid[m.N];length=50)
m.w_grid
m.a_grid
m.markovprocess[1]


xs = 1:0.2:5
f(x) = log(x)
A = [f(x) for x in xs]
interp_linear = LinearInterpolation(xs, A)
interp_linear(3) # exactly log(3)
interp_linear(3.1) # approximately log(3.1)
interp_linear(0.9) # outside grid: error
interp_linear_extrap = LinearInterpolation(xs, A,extrapolation_bc=Line()) 
interp_linear_extrap(0.9) # outside grid: linear extrapolation


filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/dolo/examples/models/consumption_savings_iid.yaml"
readlines(filename)
model = yaml_import(filename)
model.symbols[:expectation]
Dolo.arbitrage(model)

dp = Dolo.discretize(model.exogenous)
dp.grid
dp.integration_nodes
dp.integration_weights
Dolo.iweight(dp,1,2) # transition probabilities
Dolo.n_inodes(dp,1) # nb of possible future states given current state i
Dolo.n_nodes(dp) # size of matrix states
Dolo.inode(dp,1,1)


filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/julia_training/consumption-savings/CS_py.yaml"
readlines(filename)
model = yaml_import(filename)
model.symbols[:states]
dp = Dolo.discretize(model.exogenous)
dp.grid
dp
dp.transitions
dp.values
Dolo.get_integration_nodes(dp,1)

Dolo.iweight(dp,1,2) # transition probabilities
Dolo.n_inodes(dp,1) # nb of possible future states given current state i
Dolo.n_nodes(dp) # size of matrix states
Dolo.inode(dp,1,1) # which states j in particular, given i



filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/julia_training/consumption-savings/CS_py.yaml"
#filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/Dolo.jl/examples/models/consumption_savings_iid.yaml"
readlines(filename)
model = yaml_import(filename)


#=
funs = Dolo.direct_response_egm(model, 1, 1, 1, 1)
h = funs["expectation"]
gt = funs["half_transition"]
τ = funs["direct_response_egm"]
aτ = funs["reverse_state"]
=#

Dolo.get_equation_block(model, "arbitrage")
eqs = Dolo.get_assignment_block(model, "half_transition")
eqs

model.symbols
model.symbols[:poststates]
model.data["equations"]
Dolo.get_variables(model)
Dolo.expectation()

# liste
# SVector vector static, dont on connait la dimension
import Dolang
f = Dolo.get_factory(model, "transition")
code = Dolang.gen_generated_gufun(f)
transition = eval(code)
using StaticArrays
m_, s_, x_, p_ = model.calibration[:exogenous, :states, :controls, :parameters]
m,s,x,p = [SVector(e...) for e in model.calibration[:exogenous, :states, :controls, :parameters]]
model.symbols[:parameters][5]
model.calibration
mvec = [SVector(m...)  for i = 1:10]
mvec = cat([m' for i=1:100]...; dims=1)
transition(m,s,x,m,p)


for i in 1:length(vector_F)
    code = Dolang.gen_generated_gufun(vector_F[i])
    push!(F, eval(code))
end

for i in 1:length(vector_F)
    print(i)
end



Dolo.get_calibration(model)
Dolo.get_defined_variables(model)
Dolo.get_variables(model)
model.symbols

model.symbols[:poststates]
model.data["equations"]
Dolo.get_variables(model)

model.symbols
model.symbols[:poststates]
model.data["equations"]
Dolo.get_variables(model)



mvec = [SVector(m...)  for i = 1:10]
mvec = cat([m' for i=1:10]...; dims=1)
transition(m,s,x,m,p)

# testing n_inodes for different processes

# IID 
mutable struct DiscretizedIIDProcess <: AbstractDiscretizedProcess
    # nodes::Matrix{Float64}
    grid::EmptyGrid
    integration_nodes::Matrix{Float64}
    integration_weights::Vector{Float64}
end

DiscretizedIIDProcess(x, w) = DiscretizedIIDProcess(EmptyGrid{size(x,2)}(), x, w)

n_nodes(dp::DiscretizedIIDProcess) = 0
n_inodes(dp::DiscretizedIIDProcess, i::Int) = size(dp.integration_nodes, 1)
inode(dp::DiscretizedIIDProcess, i::Int, j::Int) = dp.integration_nodes[j, :]
iweight(dp::DiscretizedIIDProcess, i::Int, j::Int) = dp.integration_weights[j]
node(dip::DiscretizedIIDProcess, i::Int) = fill(NaN, n_inodes(dip, 1))

# MC 
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/Dolo.jl/examples/models/rbc_mc.yaml"
readlines(filename)
model = yaml_import(filename)

shock = Dolo.get_exogenous(model)
dp = Dolo.discretize(shock)
typeof(dp)
dp.values
dp.transitions
dp.grid
dp.grid.nodes # nodes only, no n
dp.grid.nodes[1]

Dolo.n_nodes(dp) # size of matrix states
Dolo.n_inodes(dp,1)
Dolo.iweight(dp,1,1)
Dolo.default_index(dp)
dp
Dolo.inode(dp,1,1)
Dolo.node(dp,1)

# AR
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/Dolo.jl/examples/models/neoclassical.yaml"
readlines(filename)
model = yaml_import(filename)

shock = Dolo.get_exogenous(model)
dp = Dolo.discretize(shock)
typeof(dp)
dp.values
dp.transitions
dp.grid
dp.grid.n # exog grid 
dp.grid.nodes
dp.grid.nodes[1]
dp.integration_nodes
dp.integration_nodes[5]
dp.integration_weights[1]

Dolo.n_nodes(dp) # size of matrix states
Dolo.n_inodes(dp,1) # given in node 1 you can go to # nodes
Dolo.iweight(dp,1,1)
Dolo.default_index(dp)
dp
Dolo.inode(dp,1,1)
Dolo.node(dp,1)


### integration_nodes for IID

exogrid = discretize[1].exo
exogrid.n

endogrid = discretize[1].endo
endogrid.n
grid_endo = range(grid_endo.min, grid_endo.max; length = grid_endo.n)
grid_endo.min
typeof(grid_endo)
Dolo.n_nodes(grid_endo)

N_m = max(Dolo.n_nodes(grid_exo),1) 


grid, dprocess = Dolo.discretize(model)
dprocess
grid_endo = grid.endo
grid_exo = grid.exo

N_m = max(Dolo.n_nodes(grid_exo),1) 
Dolo.inode(dprocess,1,1)
grid.endo.n
grid.endo.nodes
Dolo.nodes(grid_endo)
grid_fixed = grid_endo

nodes = Dolo.inode(dprocess,1,7)
Dolo.nodes(grid_exo)

####

sigma = Array{Float64}(undef, 1, 1)
    sigma[1,1] = 0.01
    ar1 = Dolo.VAR1(ρ = 0.95, Σ = sigma)
    mc = Dolo.discretize(ar1; n=9)
    states = mc.values
    transitions = mc.transitions

    ###

    shock = Dolo.get_exogenous(model)
    Dolo.get_domain(shock) # does not work
    grid, dprocess = Dolo.discretize(model)
Dolo.get_integration_nodes(dprocess) # does not work 
dprocess.integration_nodes #works!
Dolo.inode(dprocess,1,1) #pb is it only takes one value and not the vector of integration nodes
x = [Dolo.inode(dprocess, 1,i) for i=1:size(dprocess.integration_nodes,1)] # actually useless, line below easier
dprocess.integration_nodes
dprocess.integration_weights
Dolo.iweight(dp,1,1)

function discretize(mvn::MvNormal; n=5::Union{Int, Vector{Int}}) #does not work
    x, w = QE.qnwnorm(n, mvn.μ, mvn.Σ)
    DiscretizedIIDProcess(x, w)
end


p = SVector(model.calibration[:parameters]...)
x = SVector(model.calibration[:controls]...)

# previously in egmDolo.jl

arb = Dolo.get_equation_block(model, "arbitrage")
tran = Dolo.get_assignment_block(model, "transition")
h = Dolo.get_assignment_block(model, "expectation")
g = Dolo.get_assignment_block(model, "half_transition")
τ = Dolo.get_assignment_block(model, "direct_response_egm")
aτ = Dolo.get_assignment_block(model, "reverse_state")

dp = Dolo.discretize(shock, n = grid[:exo][:n])
typeof(dp)
dp.grid
dp.integration_nodes
dp.integration_weights
Dolo.iweight(dp,1,1) # transition probabilities
Dolo.n_inodes(dp,1) # nb of possible future states given current state i
size_states = Dolo.n_nodes(dp) # size of matrix states


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



# wrong arguments, get_factory, gen_gufun, RECIPES #

## TO RUN BEFORE THE REST 
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/Dolo.jl/examples/models/consumption_savings_iid.yaml"
readlines(filename)
model = yaml_import(filename)
# recursively make all keys at any layer of nesting a symbol
# included here instead of util.jl so we can call it on RECIPES below
_symbol_dict(x) = x
_symbol_dict(d::AbstractDict) =
    Dict{Symbol,Any}([(Symbol(k), _symbol_dict(v)) for (k, v) in d])
const src_path = dirname(@__FILE__)
const pkg_path = dirname(src_path)
Pkg.add("YAML")
import YAML; using YAML: load_file, load
# add recipes.yaml wherever the next line indiccates the error 
const RECIPES = _symbol_dict(load_file(joinpath(src_path, "recipes.yaml")))
specs = RECIPES[:dtcc][:specs]

# 6.07 fail occurs here #
# I run only the part that deals with argument in get_factory()
Symbol("half_transition") in keys(specs) # contradiction here, works here but not in egmDolo
Dolo.get_factory(model, "half_transition") # at least sign that the good file is being used 
recipe = specs[Symbol("half_transition")]
symbols = Dolo.get_symbols(model)
arguments = Dolo.OrderedDict(
            Symbol(l[3]) => [Dolang.stringify(e,l[2]) for e in symbols[Symbol(l[1])]]
            for l in recipe[:eqs] if !(l[1]=="parameters")
        )
# No problem, so it must be above in the function, or I am calling the function from the wrong file 

# If I run the entire function get_factory(model::Model, eq_type::String) instead, I get:
F_g = Dolo.get_factory(model, "half_transition")
# moving specs = RECIPES[:dtcc][:specs] outside if-loop
F_g = Dolo.get_factory(model, "half_transition")
# still... 
code_g = Dolang.gen_generated_gufun(F_g)




# interpolating #

struct MyDR 
    itp::Array{Any}
end

mydrtest = MyDR([1,2])
mydrtest.itp
(mydr::MyDR)(i,s) = MyDR.itp[i](s)
mydrtest.itp[1](s)

(mydr::MyDR)(1,1)

LinearInterpolation(s0, s0; extrapolation_bc=Line())
(mydr::MyDR)(1,1) = MyDR(LinearInterpolation(s0, s0; extrapolation_bc=Line()))
(mydr::MyDR)(i,s) = mydr.itp[i](s)


x = [1.0, 2.0, 4.0]
y = [[1.0,2.0,3.0]  [4.0,5.0,6.0]] 
itptest = Vector{Any}(undef,2)
itptest[1] = LinearInterpolation(x, y[:,1], extrapolation_bc = Line())
itptest[2] = LinearInterpolation(x, y[:,2], extrapolation_bc = Line())
res = MyDR(itptest)
mydr(i,s) = res.itp[i](s)
typeof(mydr)
mydr(1,x)
mydr(2,x)

# now function or vector ?




struct MyDR 
    itp::Vector{Any} # length i 
end

MyDR(LinearInterpolation(x, y[:,1], extrapolation_bc = Line())



(mydr::MyDR)(i,s) = mydr.itp[i](s)
mydr.itp[i]


for i in 1:2
    itptest[i] = LinearInterpolation(x, y[:,i], extrapolation_bc = Line())
    itptest[i]
end 

# end of file egmDolo

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

SVector(1,2,3...)
SVector{2,Any}(x1,x2,...)
consumption_a(φ)   #Vector{SVector}#

function φ(i, s::Vector{T})::SVector{T} where T
    # s: [w]
    # x: [c]
    
end

