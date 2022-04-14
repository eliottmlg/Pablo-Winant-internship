using Dolo
filename = "consumption_savings.jl"
readlines(filename)
model = yaml_import(filename)
residuals(model)
