using Dolo
using Pkg
Pkg.add("AxisArrays")
Pkg.add("SimplePlots")
Pkg.add("Plots")
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
tab = Dolo.tabulate(model, sol.trace.trace[4], :k)
sol.dr
dr = sol.trace.trace

# when shock is iid, shock is not a state variable, vector of states is shorter, need different indexing
plt = Plots.plot()
for i=1:length(sol.trace.trace)
    tab = Dolo.tabulate(model, sol.trace.trace[i], :k)
    Plots.plot!(plt, tab[:k], tab[:i]; legend=false)
end

plt
#plotting IRF of capital after an exogenous productivity shock
IRF = Dolo.response(model, sol.dr, :e_z)
plt = Plots.plot()
Plots.plot!(plt, 1:40, IRF[:i], label = "IRF of capital")
Plots.plot!(plt, legend = :topright)



tab = tabulate(:k, :i)
Plot(sol.dr.[:i], sol[:k])
