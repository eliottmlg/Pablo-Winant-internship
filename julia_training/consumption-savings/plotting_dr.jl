using Dolo
using Pkg
Pkg.add("AxisArrays")
Pkg.add("SimplePlots")
using SimplePlots
using AxisArrays
using Plots

# declaring the model used  
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/dolo/examples/models/rbc_iid.yaml"
readlines(filename)
model = yaml_import(filename)

# computing the solution 
sol = Dolo.time_iteration(model; trace = true)

# plotting the decision rule (investment) as a function of the state variable (capital)
tab = Dolo.tabulate(model, sol.dr, :k)

#plotting IRF of capital after an exogenous productivity shock
IRF = Dolo.response(model, sol.dr, :e_z)
plt = Plots.plot()
Plots.plot!(plt, 1:40, IRF[:k], label = "IRF of capital")
Plots.plot!(plt, legend = :topright)



tab = tabulate(:k, :i)
Plot(sol.dr.[:i], sol[:k])
