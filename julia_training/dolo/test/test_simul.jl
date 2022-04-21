

# Pkg.build("QuantEcon")
import Dolo
path = path = Dolo.pkg_path

using AxisArrays



filename = joinpath(path,"examples","models","rbc_dtcc_iid_ar1.yaml")
model = Dolo.yaml_import(filename)
@time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)

driving_process = Dolo.simulate(model2.exogenous, 10, 40, m0)
sim1 = Dolo.simulate(model, dr, driving_process)

sim = Dolo.simulate(model, dr; N= 1000, T=50)
# Ac= cat(model.symbols[:exogenous], model.symbols[:states], model.symbols[:controls]; dims=1)
# ll=[Symbol(i) for i in Ac]
# AA= AxisArray(sim, Axis{:N}(1:1000), Axis{:V}(ll), Axis{:T}(1:50))

Tvalues = linspace(1, sim[Axis{:T}][end], sim[Axis{:T}][end])

a_sort = sort(collect(sim[Axis{:V}(:k)]),1)
a_median = a_sort[Int(size(sim,1)*0.5),:]


import PyPlot
plt = PyPlot;

for i in 1:sim[Axis{:N}][end]
    plt.plot(Tvalues, sim[Axis{:N}(i), Axis{:V}(:k)], color="red", alpha = .05)
end
plt.title("Simulations");
plt.xlabel("Horizon");
plt.ylabel("Capital");


plt.plot(Tvalues, a_median, marker="o")




#########################################################################

##################  IRF         ###########################



s0 = model.calibration[:states]
e0 = model.calibration[:exogenous]
irf = Dolo.response(model, dr, s0, :e_z)
Dolo.response(model, dr, s0, [0.3])


############## 2 shocks
filename2 = joinpath(path,"examples","models","rbc_dtcc_iid_2ar1.yaml")
model2 = Dolo.yaml_import(filename2)
@time dr2 = Dolo.time_iteration(model2, verbose=true, maxit=10000)
s0 = model2.calibration[:states]
e0 = model2.calibration[:exogenous]

irf=Dolo.response(model2, dr2, s0, [0, 0.3]; T=40)

irf[:k]
irf2=Dolo.response(model2, dr2, s0, :e_z)

Dolo.response(model2, dr2, s0, :e_d)


#########################################################################

# filename = joinpath(path,"examples","models","rbc.yaml")
# model = Dolo.yaml_import(filename)
# # model2 = Dolo.yaml_import(filename2)
# N = 1
# T=40
# @time dr = Dolo.time_iteration(model, verbose=true, maxit=10000)
# s0 = model.calibration[:states]
# m0 = model.calibration[:exogenous]
# Dolo.response(model, dr, s0, [0.2])
# Dolo.response(model, dr, s0, [0.4])




############
# Ploting


# irf = Dolo.response(model, dr, e0; horizon = 100)


index = findfirst(isequal(:k), model.symbols[:states])

kirf = irf[:,3]
iirf = irf[:,2]
nirf = irf[:,4]
zirf = irf[:,5]
horizon=100
time = linspace(0,horizon-1,horizon)
using Gadfly

plot(x=time, y=kirf, Geom.point, Geom.line,
     Guide.xlabel("horizon"), Guide.ylabel("Capital"), Guide.title("IRF"))
plot(x=time, y=nirf, Geom.point, Geom.line,Guide.xlabel("horizon"),
     Guide.ylabel("Hours"), Guide.title("IRF"))
plot(x=time, y=iirf, Geom.point, Geom.line, Guide.xlabel("horizon"),
   Guide.ylabel("Investments"), Guide.title("IRF"))
plot(x=time, y=zirf, Geom.point, Geom.line, Guide.xlabel("horizon"),
      Guide.ylabel("AR1"), Guide.title("IRF"))
















# test simulation
res = Dolo.simulation(model, sigma, dr,s0, n_exp, horizon, seed, zeros(0, 0))
res_long = Dolo.simulation(model, sigma, dr,s0, n_exp=0, horizon=1000, seed=42)

res = Dolo.simulation(model, sigma, dr,s0)
res = Dolo.simulation(model, sigma)

kvec = res_long[:,2,:]
ivec = res_long[:,4,:]
nvec = res_long[:,3,:]
horizon=1000
time = linspace(0,horizon-1,horizon)
using Gadfly

plot(x=time, y=kvec, Geom.point, Geom.line,
     Guide.xlabel("horizon"), Guide.ylabel("Capital"), Guide.title("Simulations"))
plot(x=time, y=nvec, Geom.point, Geom.line,Guide.xlabel("horizon"), Guide.ylabel("Hours"), Guide.title("Simulations"))
plot(x=time, y=ivec, Geom.point, Geom.line, Guide.xlabel("horizon"), Guide.ylabel("Investments"), Guide.title("Simulations"))


O = Theme(
    default_color=colorant"orange")

plot(layer(x=time, y=nvec, Geom.point, Geom.line, O),
     layer(x=time, y=ivec, Geom.point, Geom.line),
     Guide.xlabel("horizon"),
     Guide.ylabel("Hours (orange), Investment"),
     Guide.title("Simulations"))


# verbose=true
# verbose && @printf "%-8s%-10s%-10s%-10s%-5s\n" "t" model.symbols[:states][1] model.symbols[:states][2] model.symbols[:controls][1] model.symbols[:controls][2]
# verbose && println(repeat("-", 35))


#verbose && @printf "%-8s%-10s%-10s%-10s%-5s\n"  t round(s[1],2) round(s[2],2) round(x[1],2) round(x[2],2)
