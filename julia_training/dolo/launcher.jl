
# using Pkg
#Pkg.add("Optim")
using Dolo

filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/julia_training/dolo/consumption_savings.yaml"
readlines(filename)
model = yaml_import(filename)
residuals(model)
Base.show(model)
dr_pert = perturb(model)
dr_global = time_iteration(model)
    

