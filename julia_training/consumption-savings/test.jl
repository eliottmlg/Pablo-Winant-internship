
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

