using Dolo
filename = "../dolo/consumption_savings.yaml"
readlines(filename)
model = yaml_import(filename)
residuals(model)
