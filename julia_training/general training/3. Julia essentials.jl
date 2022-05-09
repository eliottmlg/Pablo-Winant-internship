## declaring julia file 
# add .jl

# new terminal, going from normal terminal to julia terminal 
# ctrl+shift+p opens palette 
# julia: start repl starts julia
# exit repl with backspace + exit()

using LinearAlgebra, Statistics

# integers and floats
# @show, im, typeof(), $

# strings
# split(s), replace(s, "surf"=>"ski"), strip, 
match(r"(\d+)", "Top 10")  # find digits in string

# Containers
x = ("foo", "bar")
x = ("foo", "bar")
y = ("foo", 2)
typeof(x), typeof(y)

x = ("foo", 1,)
y = ("foo")
typeof(x), typeof(y)

# common data types

# 3.6 broadcasting
