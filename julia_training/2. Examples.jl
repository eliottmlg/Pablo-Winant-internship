# https://julia.quantecon.org/getting_started_julia/julia_by_example.html

# exercise 2.4.1
# building factorial function 

function factorial2(n)
    initfactorial=1
    factorial_old=initfactorial
    for i in 1:n
       factorial_new=i*factorial_old
       factorial_old=factorial_new
    end
    return(value=factorial_old)
end 
sol = factorial2(20) # computing 20!, max factorial that can be computed with Julia is 20

#solution 
function factorial2(n)
    k = 1
    for i in 1:n
        k *= i  # or k = k * i
    end
    return k
end

factorial2(4)

# exercise 2 programming a binomial variable
function binomial_rv(n, p)
    y = Vector{Float64}(undef, n)
    Y = 0
    for i in 1:n 
        y[i] = rand(1)[1] 
        if y[i] >= p
            Y = Y + 1
        end
    end
    return Y
end
binomial_rv(10,0.50)

# exercise 3 approximation to pi using monte-carlo
n = 100000
radius = Vector{Float64}(undef, n)
inB = 0.0
for i in 1:n
  U = rand(2)
  radius[i] = (U[1]^2 + U[2]^2)^(1/2)
  if radius[i] <= 1
    inB = inB + 1.0
  elseif radius[i] > 1
    inB = inB
  end
end
ProbB = inB/n # probB = area of B 
area_circle = ProbB*4.0
pi_calculated = area_circle

# exercise 7 finding first passage time of random walk

t_ = Vector{Float64}(undef, 100)
t_0 = 1
x = Vector{Float64}(undef, 201)
x[1] = 0
x[201] = 0
e = Vector{Float64}(undef, 199)
alpha = 1.0
sigma = 0.2

for j in 1:100
for i in 1:199
    e[i] = randn(1)[1]
end
for i in 1:199 
    x[i+1] = alpha*x[i] + sigma*e[i]
    if x[i+1] <= 0 
        t_0 = i+1.0
        break
    end
end
t_[j] = t_0
end

using Plots
index = 1:200
plot(index, x)



jean baptiste michau
pierre edouard collignon
savoir resoudre model consumption saving (revenue fluctue chain markov etat haut etat bas, epargner a certain taux dinteret, max la sum des consommation avec les gamma)
ABOSULEMENT: value function iteration
simple, un shock revenue exogene 
DANS UNE SEMAINE,
cours demain matin bachelor X, 8h a midi, introduction a julia, modele epidemiology, slides sur la convergence des series, partie programmation utile
COURS: site github 
mie37 blob  master notebook 

