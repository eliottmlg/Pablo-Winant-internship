
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


exogenous = UNormal(;σ=0.1)
dp = discretize(exogenous)
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
Dolo.iweight(dp,1,1) # transition probabilities
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

Dolo.get_equation_block(model, "half_transition")


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

transition(m,s,x,m,p)

m_, s_, x_, p_ = model.calibration[:exogenous, :states, :controls, :parameters]
m,s,x,p = [SVector(e...) for e in model.calibration[:exogenous, :states, :controls, :parameters]]


mvec = [SVector(m...)  for i = 1:10]
mvec = cat([m' for i=1:100]...; dims=1)

transition(m,s,x,m,p)

