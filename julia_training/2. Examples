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

# exercise 2 
function binomial_rv(n, p)
    
    
end



