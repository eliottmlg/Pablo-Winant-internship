
# Setup
using Pkg
using Dolo 

# Consumption savings Dolo
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/julia_training/consumption-savings/consumption_savings.yaml"
# Consumption savings iid Dolo
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/dolo/examples/models/consumption_savings_iid.yaml"
readlines(filename)
model = yaml_import(filename)
residuals(model)
Base.show(model)
dr_pert = perturb(model)
dr_global = Dolo.time_iteration(model) 
dr_global

dr_global_ITI = Dolo.improved_time_iteration(model)



