using Dolo
using Pkg
Pkg.add("AxisArrays")
Pkg.add("SimplePlots")
Pkg.add("Plots")
using SimplePlots
using AxisArrays
using Plots
using LaTeXStrings

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

# Consumption-savings model with iid shock
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/dolo/examples/models/consumption_savings_iid.yaml"
readlines(filename)
model = yaml_import(filename)
sol = Dolo.time_iteration(model; trace = true)
# sol.dr 
sol.dprocess
decisionrule = sol.trace.trace

xvec = range(0,10;length=100)
dr = decisionrule
tab = Dolo.tabulate(model, dr[1], :w)
w_grid = tab[:w]

function convergenceDR(decisionrule)
    dr = decisionrule
    tab = Dolo.tabulate(model, dr[1], :w)
    w_grid = tab[:w]
    plt = Plots.plot()
    plot!(plt, w_grid, w_grid; ylims = (0,2))
    plot!(plt,  w_grid, tab[:c], color = RGBA(0,0,0,1), label = L"initial condition $c(w) = constant$")
    for i=2:10:length(dr)
        tab = Dolo.tabulate(model, dr[i], :w)
        Plots.plot!(plt, tab[:w], tab[:c]; label = [i], legend=false)
    end
    plot!(plt, xlabel = "Wealth", ylabel = "Consumption")
    plot!(plt, legend = :bottomright)
end
convergenceDR(decisionrule)


