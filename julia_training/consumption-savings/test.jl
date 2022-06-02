
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