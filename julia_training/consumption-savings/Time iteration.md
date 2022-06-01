(to open preview window: ctrl+shift+v)

# Time iteration algorithm

--- 

## 1.Useful markdown commands 
1. item 1 
2. item 2
3. item 3 

*this test is italic* and **bold** 

yes 

[link](https://www.youtube.com/watch?v=DLLrcr9u_XI)

```
to add code 
```

1st Header|2nd Header|3rd Header
---|:---:|---: 
col 1 is|left-aligned|1
col 2 is|center-aligned|2
col 3 is|right-aligned|3

$\eta < \tau_{\eta}$

## 2.Description of the TI algorithm
The TI algorithm consists in applying the same function (the Coleman operator) repeatedly to the policy function (in general consumption, which is expressed in terms of the states) until the equilibrium condition of the model (in general the Euler equation) is met.

Thus the necessary and sufficient criterion for a policy function to be a solution of the model is that it is the policy function such that $\epsilon<\tau_\epsilon$. This means that the difference between the LHS and the RHS of the Euler equation is approximately zero.

One can add a secondary criterion for determining whether the algorithm has converged or not, namely if $\eta<\tau_\eta$. This means that the difference between the policy function computed in the previous iteration and that computed in the current iteration is approximately zero.

We will be testing the TI algorithm on the consumption_savings_iid.yaml file.

## 3.Tracking changes to make the TI output more user-friendly
1. We start with a TimeIterationResult object of the form:
```r
   mutable struct TimeIterationResult
    dr::AbstractDecisionRule
    iterations::Int
    complementarities::Bool
    dprocess::AbstractDiscretizedProcess
    x_converged::Bool
    x_tol::Float64  
    err::Float64
    log::IterationLog
    trace::Union{Nothing,IterationTrace}
    end
```

and a function converged() of the form:

```r
converged(r::TimeIterationResult) = r.x_converged
```
The TI algorithm would yield the following output:

```r
function Base.show(io::IO, r::TimeIterationResult)
    @printf io "Results of Time Iteration Algorithm\n"
    @printf io " * Complementarities: %s\n" string(r.complementarities)
    @printf io " * Discretized Process type: %s\n" string(typeof(r.dprocess))
    @printf io " * Decision Rule type: %s\n" string(typeof(r.dr))
    @printf io " * Number of iterations: %s\n" string(r.iterations)
    @printf io " * Convergence: %s\n" converged(r)
    @printf io "   * |x - x'| < %.1e: %s\n" r.x_tol r.x_converged
end
```

and the following TI algorithm:

```r
function time_iteration(model;
    dr0=Dolo.ConstantDecisionRule(model),
    discretization=Dict(),
    interpolation=:cubic,
    verbose=true,
    details=true,
    ignore_constraints=false,
    trace = false,
    tol_η = 1e-8,
    tol_ε = 1e-8,
    maxit=1000
)
    
    F = Euler(model; discretization=discretization, interpolation=interpolation, dr0=dr0,  ignore_constraints=ignore_constraints)

    complementarities = false


    if interpolation != :cubic
        error("Interpolation option ($interpolation) is currently not recognized.")
    end

    ti_trace = trace ? IterationTrace([x0]) : nothing


    z0 = deepcopy(F.x0)
    z1 = deepcopy(z0)

    log = IterationLog(
        it = ("n", Int),
        err =  ("ϵₙ=|F(xₙ,xₙ)|", Float64),
        sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
        lam = ("λₙ=ηₙ/ηₙ₋₁",Float64),
        elapsed = ("Time", Float64)
    )

    initialize(log, verbose=verbose; message="Time Iteration")

    local err_η, z0, z1

    err_η_0 = NaN

    it = 0

    while it<=maxit

        it += 1

        t1 = time_ns()

        r0 = F(z0, z0; set_future=true)

        err_ε= norm(r0)
        if err_ε<tol_ε
            break
        end

        fun = u->F(u, z0; set_future=false)
        dfun = u->df_A(F,u, z0; set_future=false)

        sol = newton2(
            fun, dfun,
            z0, verbose=false
        )

        z1 = sol
        δ = z1 - z0

        trace && push!(ti_trace.trace, z1)

        err_η = norm(δ)
        gain = err_η / err_η_0
        err_η_0 = err_η

        if err_η<tol_η
            break 
        end

        # z0.data[:] .= z1.data
        z0 = z1

        elapsed = time_ns() - t1

        elapsed /= 1000000000

        append!(log; verbose=verbose, it=it, err=err_ε, sa=err_η, lam=gain, elapsed=elapsed)

    end

    finalize(log, verbose=verbose)


    res = TimeIterationResult(F.dr.dr, it, complementarities, F.dprocess, err_η<tol_η, tol_η, err_η, log, ti_trace)
    return res

end
```
the output is the following

```
108 |   1.3609e-08 |   1.5677e-08 |   8.0463e-01 |   5.2709e-02
     109 |   1.0950e-08 |   1.2615e-08 |   8.0465e-01 |   5.4535e-02
------------------------------------------------------------------
Results of Time Iteration Algorithm
 * Complementarities: false
 * Discretized Process type: Dolo.DiscretizedIIDProcess
 * Decision Rule type: Dolo.CubicDR{Dolo.EmptyGrid{1}, Dolo.UCGrid{1}, 1, 1}
 * Number of iterations: 110
 * Convergence: false
   * |x - x'| < 1.0e-08: false
```

according to this algorithm, convergence is achieved when the distance between two policy function computed is zero. Here, "Convergence: false" even though the algorithm stopped. The reason is that   what is displayed is not the current iteration (110th) but the precedent.

We modify the function time_iteration(). 

```r

        append!(log; verbose=verbose, it=it, err=err_ε, sa=err_η, lam=gain, elapsed=elapsed)

        if err_ε<tol_ε       # change
            break
        end

        if err_η<tol_η       # change
            break 
        end

    end

    finalize(log, verbose=verbose)


    res = TimeIterationResult(F.dr.dr, it, complementarities, F.dprocess, err_η<tol_η, tol_η, err_η, log, ti_trace)
    return res

end
```

which yields the following:

```
109 |   1.0950e-08 |   1.2615e-08 |   8.0465e-01 |   6.4835e-02
     110 |   8.8111e-09 |   1.0151e-08 |   8.0466e-01 |   5.5368e-02
------------------------------------------------------------------
Results of Time Iteration Algorithm
 * Complementarities: false
 * Discretized Process type: Dolo.DiscretizedIIDProcess
 * Decision Rule type: Dolo.CubicDR{Dolo.EmptyGrid{1}, Dolo.UCGrid{1}, 1, 1}
 * Number of iterations: 110
 * Convergence: false
   * |x - x'| < 1.0e-08: false
```

so we fixed the problem of the wrong iteration being displayed. 

2. We change the criterion for convergence, which we now set as whether the Euler equation is satisfied. 

we first remove x_converged, which indicates if the distance between two policy function is close to zero, and can be misinterpreted with the event the Euler equation is satisfied.

```r
  mutable struct TimeIterationResult
    dr::AbstractDecisionRule
    iterations::Int
    complementarities::Bool
    dprocess::AbstractDiscretizedProcess
    ϵ::Float64 # distance between RHS and LFS of Euler equation, previously err_ε
    τ_ϵ::Float64 # tolerance for ϵ, previously tol_ε
    η::Float64 # distance between two policy functions computed, previously err or err_η
    τ_η::Float64 # tolerance for η, previously x_tol or tol_η
    log::IterationLog
    trace::Union{Nothing,IterationTrace}
end
```

we then redefine the converged() function so that is it TRUE when the TI algorithm finds a solution, that is when the Euler equation is met. For this, need the following:
```r
converged(r::TimeIterationResult) = r.ϵ<r.τ_ϵ 
```
The TI algorithm would yield the following output:

```r
function Base.show(io::IO, r::TimeIterationResult)
    @printf io "Results of Time Iteration Algorithm\n"
    @printf io " * Complementarities: %s\n" string(r.complementarities)
    @printf io " * Discretized Process type: %s\n" string(typeof(r.dprocess))
    @printf io " * Decision Rule type: %s\n" string(typeof(r.dr))
    @printf io " * Number of iterations: %s\n" string(r.iterations)
    @printf io " * Convergence: %s\n" converged(r) 
end
```
where we remove |x - x'| for redundancy.

with the following algorithm

```r
function time_iteration(model;
    dr0=Dolo.ConstantDecisionRule(model),
    discretization=Dict(),
    interpolation=:cubic,
    verbose=true,
    details=true,
    ignore_constraints=false,
    trace = false,
    τ_ϵ = 1e-8,
    τ_η = 1e-8,
    maxit=500
)
    
    F = Euler(model; discretization=discretization, interpolation=interpolation, dr0=dr0,  ignore_constraints=ignore_constraints)

    complementarities = false


    if interpolation != :cubic
        error("Interpolation option ($interpolation) is currently not recognized.")
    end

    ti_trace = trace ? IterationTrace([x0]) : nothing


    z0 = deepcopy(F.x0)
    z1 = deepcopy(z0)

    log = IterationLog(
        it = ("n", Int),
        err =  ("ϵₙ=|F(xₙ,xₙ)|", Float64),
        sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
        lam = ("λₙ=ηₙ/ηₙ₋₁",Float64),
        elapsed = ("Time", Float64)    
    )                                 

    initialize(log, verbose=verbose; message="Time Iteration")

    local η, z0, z1

    η_0 = NaN

    it = 0

    while it<=maxit

        it += 1

        t1 = time_ns()

        r0 = F(z0, z0; set_future=true)

        ϵ = norm(r0)

        fun = u->F(u, z0; set_future=false)
        dfun = u->df_A(F,u, z0; set_future=false)

        sol = newton2(
            fun, dfun,
            z0, verbose=false
        )

        z1 = sol
        δ = z1 - z0

        trace && push!(ti_trace.trace, z1)

        η = norm(δ)
        gain = η / η_0
        η_0 = η


        # z0.data[:] .= z1.data
        z0 = z1

        elapsed = time_ns() - t1

        elapsed /= 1000000000

        append!(log; verbose=verbose, it=it, err=ϵ, sa=η, lam=gain, elapsed=elapsed)

        
        if ϵ<τ_ϵ      
           break                 
       end  

        if η<τ_η     
            break 
         end

    end

    finalize(log, verbose=verbose)


    res = TimeIterationResult(F.dr.dr, it, complementarities, F.dprocess, ϵ, τ_ϵ, η, τ_η, log, ti_trace)  # line 313
  return res

end
```
however it encounters an error

```r
  110 |   8.8111e-09 |   1.0151e-08 |   8.0466e-01 |   5.8307e-02
------------------------------------------------------------------
ERROR: UndefVarError: err not defined
Stacktrace:
 [1] time_iteration(model::Model{Symbol("##292"), Dolo.ProductDomain{2, Dolo.EmptyDomain{1}, Dolo.CartesianDomain{1}}}; dr0::Dolo.ConstantDecisionRule{1}, discretization::Dict{Any, Any}, interpolation::Symbol, verbose::Bool, details::Bool, ignore_constraints::Bool, trace::Bool, τ_ϵ::Float64, τ_η::Float64, maxit::Int64)    
   @ Dolo c:\Users\t480\GitHub\Pablo-Winant-internship\Dolo\src\algos\time_iteration.jl:313
```

It seems it does not make sense to include the tolerance of epsilon and eta in the structure TimeIterationResult, instead it would make sense to include a convergence criterion, for both epsilon and eta. 

Still the same problem, variable undefined. In Julia, a variable x requires a new binding each time a loop is executed. The was a problem with the scope of the variable err_$\epsilon$, which was only defined within the ``while`` function, and was thus not passed onto the functions out of the loop, such as ``res``. 

We solve the problem by defining the variable err_$\epsilon$ outside the ``while`` function, using ``local err_$\epsilon$, err_$\eta$``. We make sure to use the good font for greek letters. 

```r
    local err_ε, err_η, z0, z1 #

    log = IterationLog(
        it = ("n", Int),
        err =  ("εₙ=|F(xₙ,xₙ)|", Float64),
        sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
        lam = ("λₙ=ηₙ/ηₙ₋₁",Float64),
        elapsed = ("Time", Float64)
    )

    initialize(log, verbose=verbose; message="Time Iteration")

    err_η_0 = NaN

    it = 0

    while it<=maxit

        it += 1

        t1 = time_ns()

        r0 = F(z0, z0; set_future=true)

        err_ε=norm(r0)

        fun = u->F(u, z0; set_future=false)
        dfun = u->df_A(F,u, z0; set_future=false)

        sol = newton2(
            fun, dfun,
            z0, verbose=false
        )

        z1 = sol
        δ = z1 - z0

        #trace && push!(ti_trace.trace, z1)

        err_η = norm(δ)
        gain = err_η / err_η_0
        err_η_0 = err_η

        # z0.data[:] .= z1.data
        z0 = z1

        elapsed = time_ns() - t1

        elapsed /= 1000000000

     append!(log; 
        verbose=verbose, 
        it=it, 
        err=err_ε, 
        sa=err_η, 
        lam=gain, 
        elapsed=elapsed
     )

     if err_ε<tol_ε     
         break
     end

     if err_η<tol_η      
         break 
     end

end

 finalize(log, verbose=verbose)


 res = TimeIterationResult(
     F.dr.dr, 
     it, 
     complementarities, 
     F.dprocess, 
     err_ε, 
     err_η, 
     tol_ε, 
     tol_η, 
     log, 
     ti_trace
     )

 return res

end
```
3. We modify the code to display the minimal number of iterations left before eta meets its tolerance level.
   We replace all ``log`` with ``Log`` in the file "time_iteration.jl" in order to use the logarithm function to calculate $p$. 

   where
   $p \ge \frac{log(\tau_\eta/\eta_n)}{log(\lambda)}$. 

We also modify the code to show the time remaining, before eta meets its tolerance level. 

$TimeLeft = \frac{time_n + time_{n-1} + ... + time_{n-4}}{5} * p$

```r

mutable struct TimeIterationResult
    dr::AbstractDecisionRule
    iterations::Int
    complementarities::Bool
    dprocess::AbstractDiscretizedProcess
    ϵ::Float64 
    η::Float64
    τ_ϵ::Float64 
    τ_η::Float64
    τ_λ::Float64
    p::Int #
    Log::IterationLog
    trace::Union{Nothing,IterationTrace}
end

function time_iteration(model;
    dr0=Dolo.ConstantDecisionRule(model),
    discretization=Dict(),
    interpolation=:cubic,
    verbose=true,
    details=true,
    ignore_constraints=false,
    trace = false,
    tol_η = 1e-8,
    tol_ε = 1e-8,
    λbar = 8.0466e-01, #
    maxit=1000,
)
    
    F = Euler(model; discretization=discretization, interpolation=interpolation, dr0=dr0,  ignore_constraints=ignore_constraints)

    complementarities = false


    if interpolation != :cubic
        error("Interpolation option ($interpolation) is currently not recognized.")
    end

    ti_trace = trace ? IterationTrace([x0]) : nothing


    z0 = deepcopy(F.x0)
    z1 = deepcopy(z0)

    vector_time = rand(maxit) #

    local err_ε, err_η, z0, z1, p #

    Log = IterationLog(
        it = ("n", Int),
        err =  ("εₙ=|F(xₙ,xₙ)|", Float64),
        sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
        lam = ("λₙ=ηₙ/ηₙ₋₁", Float64),
        elapsed = ("Time", Float64),
        nb_it_before_convergence_of_x = ("Min nb of iter before ηₙ<τ_η", Int), #
        remaining_time = ("Time left before ηₙ<τ_η", Float64) #
    )

    initialize(Log, verbose=verbose; message="Time Iteration")

    err_η_0 = NaN

    it = 0

    while it<=maxit

        it += 1

        t1 = time_ns()

        r0 = F(z0, z0; set_future=true)

        err_ε=norm(r0)

        fun = u->F(u, z0; set_future=false)
        dfun = u->df_A(F,u, z0; set_future=false)

        sol = newton2(
            fun, dfun,
            z0, verbose=false
        )

        z1 = sol
        δ = z1 - z0

        #trace && push!(ti_trace.trace, z1)

        err_η = norm(δ)
        gain = err_η / err_η_0
        err_η_0 = err_η

        # z0.data[:] .= z1.data
        z0 = z1

        p = log(tol_η/err_η) / log(λbar) # 

        elapsed = time_ns() - t1

        elapsed /= 1000000000

        vector_time[it] = elapsed #

        local avg_time #

        for it in 5:maxit #
            avg_time = (vector_time[it]+vector_time[it-1]+vector_time[it-2]+vector_time[it-3]+vector_time[it-4])/5
        end

        time_left = avg_time*p # 

        p = round(Int, p) #

     append!(Log; 
        verbose=verbose, 
        it=it, 
        err=err_ε, 
        sa=err_η, 
        lam=gain, 
        elapsed=elapsed,
        nb_it_before_convergence_of_x=p, #
        remaining_time=time_left #
     )

     if err_ε<tol_ε     
         break
     end

     if err_η<tol_η       
         break 
     end

end






------------------------------------------------------------------------------------------------------------------------------
Time Iteration
------------------------------------------------------------------------------------------------------------------------------
       n | εₙ=|F(xₙ,xₙ)| | ηₙ=|xₙ-xₙ₋₁| |   λₙ=ηₙ/ηₙ₋₁ |         Time | Min nb of iter before ηₙ<τ_η | Time 
left before ηₙ<τ_η
------------------------------------------------------------------------------------------------------------------------------
       1 |   3.9481e-01 |   4.0000e-01 |          NaN |   1.3877e-01 |       81 |   3.7378e+01
       2 |   5.2508e-03 |   4.7671e-03 |   1.1918e-02 |   5.5452e-02 |       60 |   2.7919e+01
       3 |   5.2505e-03 |   4.7919e-03 |   1.0052e+00 |   1.2458e-01 |       60 |   2.7930e+01
       4 |   5.2502e-03 |   4.8172e-03 |   1.0053e+00 |   6.9226e-02 |       60 |   2.7941e+01
       5 |   5.2499e-03 |   4.8425e-03 |   1.0053e+00 |   6.4288e-02 |       60 |   2.7953e+01
       6 |   5.2495e-03 |   4.8680e-03 |   1.0053e+00 |   6.7813e-02 |       60 |   2.7964e+01
       7 |   5.2492e-03 |   4.8937e-03 |   1.0053e+00 |   1.1304e-01 |       60 |   2.7975e+01
       8 |   5.2489e-03 |   4.9195e-03 |   1.0053e+00 |   7.6218e-02 |       60 |   2.7986e+01
       9 |   5.2486e-03 |   4.9454e-03 |   1.0053e+00 |   1.4615e-01 |       60 |   2.7997e+01
      10 |   5.2483e-03 |   4.9714e-03 |   1.0053e+00 |   1.2339e-01 |       60 |   2.8009e+01
      11 |   5.2479e-03 |   4.9976e-03 |   1.0053e+00 |   2.1825e-01 |       60 |   2.8020e+01
      12 |   5.2476e-03 |   5.0240e-03 |   1.0053e+00 |   2.5142e-01 |       60 |   2.8031e+01
      13 |   5.2473e-03 |   5.0504e-03 |   1.0053e+00 |   1.8880e-01 |       60 |   2.8042e+01
      14 |   5.2470e-03 |   5.0771e-03 |   1.0053e+00 |   2.1069e-01 |       60 |   2.8054e+01
      15 |   5.2467e-03 |   5.1038e-03 |   1.0053e+00 |   1.1024e-01 |       60 |   2.8065e+01
      16 |   5.2464e-03 |   5.1307e-03 |   1.0053e+00 |   8.7310e-02 |       60 |   2.8076e+01
      17 |   5.2461e-03 |   5.1577e-03 |   1.0053e+00 |   7.3690e-02 |       61 |   2.8087e+01
      18 |   5.2458e-03 |   5.1849e-03 |   1.0053e+00 |   7.6430e-02 |       61 |   2.8098e+01
      19 |   5.2455e-03 |   5.2122e-03 |   1.0053e+00 |   8.8334e-02 |       61 |   2.8110e+01
      20 |   5.2452e-03 |   5.2397e-03 |   1.0053e+00 |   8.4405e-02 |       61 |   2.8121e+01
      21 |   5.2449e-03 |   5.2673e-03 |   1.0053e+00 |   6.7848e-02 |       61 |   2.8132e+01
      22 |   5.2446e-03 |   5.2950e-03 |   1.0053e+00 |   6.6734e-02 |       61 |   2.8143e+01
      23 |   5.2443e-03 |   5.3229e-03 |   1.0053e+00 |   9.4400e-02 |       61 |   2.8155e+01
      24 |   5.2440e-03 |   5.3510e-03 |   1.0053e+00 |   7.0300e-02 |       61 |   2.8166e+01
      25 |   5.2437e-03 |   5.3792e-03 |   1.0053e+00 |   7.2132e-02 |       61 |   2.8177e+01
      26 |   5.2434e-03 |   5.4075e-03 |   1.0053e+00 |   7.1074e-02 |       61 |   2.8188e+01
      27 |   5.2431e-03 |   5.4360e-03 |   1.0053e+00 |   1.0099e-01 |       61 |   2.8199e+01
      28 |   5.2428e-03 |   5.4646e-03 |   1.0053e+00 |   6.6415e-02 |       61 |   2.8211e+01
      29 |   5.2425e-03 |   5.4934e-03 |   1.0053e+00 |   6.4823e-02 |       61 |   2.8222e+01
      30 |   5.2422e-03 |   5.5224e-03 |   1.0053e+00 |   6.6237e-02 |       61 |   2.8233e+01
      31 |   5.2420e-03 |   5.5515e-03 |   1.0053e+00 |   9.2816e-02 |       61 |   2.8244e+01
      32 |   5.2417e-03 |   5.5807e-03 |   1.0053e+00 |   1.0534e-01 |       61 |   2.8256e+01
      33 |   5.2414e-03 |   5.6101e-03 |   1.0053e+00 |   6.4317e-02 |       61 |   2.8267e+01
      34 |   5.2411e-03 |   5.6397e-03 |   1.0053e+00 |   6.5624e-02 |       61 |   2.8278e+01
      35 |   5.2408e-03 |   5.6694e-03 |   1.0053e+00 |   9.8675e-02 |       61 |   2.8289e+01
      36 |   5.2405e-03 |   5.6993e-03 |   1.0053e+00 |   6.0696e-02 |       61 |   2.8300e+01
      37 |   5.2403e-03 |   5.7293e-03 |   1.0053e+00 |   6.8301e-02 |       61 |   2.8312e+01
      38 |   5.2400e-03 |   5.7595e-03 |   1.0053e+00 |   5.9548e-02 |       61 |   2.8323e+01
      39 |   5.2397e-03 |   5.7898e-03 |   1.0053e+00 |   6.4336e-02 |       61 |   2.8334e+01
      40 |   5.2394e-03 |   5.8203e-03 |   1.0053e+00 |   8.8881e-02 |       61 |   2.8345e+01
      41 |   5.2392e-03 |   5.8510e-03 |   1.0053e+00 |   6.0697e-02 |       61 |   2.8357e+01
      42 |   5.2389e-03 |   5.8818e-03 |   1.0053e+00 |   6.3033e-02 |       61 |   2.8368e+01
      43 |   5.2385e-03 |   5.9127e-03 |   1.0053e+00 |   6.4835e-02 |       61 |   2.8379e+01
      44 |   5.2379e-03 |   5.9435e-03 |   1.0052e+00 |   9.1779e-02 |       61 |   2.8390e+01
      45 |   5.2365e-03 |   5.9733e-03 |   1.0050e+00 |   6.5730e-02 |       61 |   2.8401e+01
      46 |   5.2323e-03 |   5.9999e-03 |   1.0044e+00 |   6.3904e-02 |       61 |   2.8410e+01
      47 |   5.2205e-03 |   6.0170e-03 |   1.0028e+00 |   7.0077e-02 |       61 |   2.8416e+01
      48 |   5.1902e-03 |   6.0109e-03 |   9.9899e-01 |   1.0120e-01 |       61 |   2.8414e+01
      49 |   5.1214e-03 |   5.9564e-03 |   9.9093e-01 |   6.0948e-02 |       61 |   2.8395e+01
      50 |   4.9847e-03 |   5.8166e-03 |   9.7652e-01 |   6.2326e-02 |       61 |   2.8344e+01
      51 |   4.7480e-03 |   5.5517e-03 |   9.5446e-01 |   6.2266e-02 |       61 |   2.8244e+01
      52 |   4.3916e-03 |   5.1380e-03 |   9.2549e-01 |   8.9154e-02 |       61 |   2.8079e+01
      53 |   3.9224e-03 |   4.5857e-03 |   8.9251e-01 |   5.7018e-02 |       60 |   2.7836e+01
      54 |   3.3768e-03 |   3.9414e-03 |   8.5950e-01 |   5.8147e-02 |       59 |   2.7513e+01
      55 |   2.8085e-03 |   3.2712e-03 |   8.2995e-01 |   6.2311e-02 |       58 |   2.7115e+01
      56 |   2.2684e-03 |   2.6364e-03 |   8.0594e-01 |   9.0571e-02 |       57 |   2.6654e+01
      57 |   1.7911e-03 |   2.0776e-03 |   7.8806e-01 |   5.9007e-02 |       56 |   2.6146e+01
      58 |   1.3920e-03 |   1.6120e-03 |   7.7587e-01 |   6.0726e-02 |       55 |   2.5604e+01
      59 |   1.0709e-03 |   1.2386e-03 |   7.6838e-01 |   6.1701e-02 |       54 |   2.5041e+01
      60 |   8.1958e-04 |   9.4696e-04 |   7.6453e-01 |   6.3140e-02 |       53 |   2.4468e+01
      61 |   6.2606e-04 |   7.2282e-04 |   7.6331e-01 |   9.0010e-02 |       51 |   2.3891e+01
      62 |   4.7851e-04 |   5.5215e-04 |   7.6389e-01 |   5.9891e-02 |       50 |   2.3316e+01
      63 |   3.6652e-04 |   4.2275e-04 |   7.6565e-01 |   6.4942e-02 |       49 |   2.2746e+01
      64 |   2.8161e-04 |   3.2472e-04 |   7.6811e-01 |   6.1935e-02 |       48 |   2.2182e+01
      65 |   2.1716e-04 |   2.5034e-04 |   7.7094e-01 |   8.4066e-02 |       47 |   2.1627e+01
      66 |   1.6809e-04 |   1.9374e-04 |   7.7391e-01 |   6.0664e-02 |       45 |   2.1080e+01
      67 |   1.3060e-04 |   1.5051e-04 |   7.7687e-01 |   5.8528e-02 |       44 |   2.0540e+01
      68 |   1.0184e-04 |   1.1736e-04 |   7.7972e-01 |   8.0356e-02 |       43 |   2.0009e+01
      69 |   7.9686e-05 |   9.1821e-05 |   7.8240e-01 |   8.7610e-02 |       42 |   1.9485e+01
      70 |   6.2546e-05 |   7.2067e-05 |   7.8487e-01 |   6.4393e-02 |       41 |   1.8968e+01
      71 |   4.9234e-05 |   5.6727e-05 |   7.8713e-01 |   7.0366e-02 |       40 |   1.8457e+01
      72 |   3.8856e-05 |   4.4767e-05 |   7.8918e-01 |   7.5161e-02 |       39 |   1.7951e+01
      73 |   3.0736e-05 |   3.5412e-05 |   7.9102e-01 |   1.0194e-01 |       38 |   1.7451e+01
      74 |   2.4364e-05 |   2.8070e-05 |   7.9266e-01 |   6.5230e-02 |       37 |   1.6954e+01
      75 |   1.9348e-05 |   2.2291e-05 |   7.9413e-01 |   6.5020e-02 |       35 |   1.6462e+01
      76 |   1.5390e-05 |   1.7731e-05 |   7.9542e-01 |   6.7362e-02 |       34 |   1.5973e+01
      77 |   1.2260e-05 |   1.4124e-05 |   7.9657e-01 |   9.2035e-02 |       33 |   1.5488e+01
      78 |   9.7781e-06 |   1.1265e-05 |   7.9758e-01 |   5.3339e-02 |       32 |   1.5005e+01
      79 |   7.8076e-06 |   8.9947e-06 |   7.9848e-01 |   6.0642e-02 |       31 |   1.4524e+01
      80 |   6.2403e-06 |   7.1891e-06 |   7.9926e-01 |   7.0158e-02 |       30 |   1.4046e+01
      81 |   4.9920e-06 |   5.7509e-06 |   7.9995e-01 |   8.0362e-02 |       29 |   1.3569e+01
      82 |   3.9964e-06 |   4.6039e-06 |   8.0056e-01 |   1.0221e-01 |       28 |   1.3094e+01
      83 |   3.2014e-06 |   3.6882e-06 |   8.0109e-01 |   6.4281e-02 |       27 |   1.2621e+01
      84 |   2.5661e-06 |   2.9562e-06 |   8.0155e-01 |   6.3921e-02 |       26 |   1.2148e+01
      85 |   2.0579e-06 |   2.3708e-06 |   8.0196e-01 |   6.3815e-02 |       25 |   1.1677e+01
      86 |   1.6511e-06 |   1.9021e-06 |   8.0231e-01 |   9.6444e-02 |       24 |   1.1207e+01
      87 |   1.3252e-06 |   1.5267e-06 |   8.0262e-01 |   6.4084e-02 |       23 |   1.0737e+01
      88 |   1.0640e-06 |   1.2258e-06 |   8.0290e-01 |   6.7109e-02 |       22 |   1.0268e+01
      89 |   8.5454e-07 |   9.8445e-07 |   8.0313e-01 |   6.3758e-02 |       21 |   9.8002e+00
      90 |   6.8649e-07 |   7.9085e-07 |   8.0334e-01 |   9.2198e-02 |       20 |   9.3326e+00
      91 |   5.5161e-07 |   6.3547e-07 |   8.0352e-01 |   6.8112e-02 |       19 |   8.8655e+00
      92 |   4.4332e-07 |   5.1071e-07 |   8.0368e-01 |   6.3140e-02 |       18 |   8.3988e+00
      93 |   3.5635e-07 |   4.1052e-07 |   8.0382e-01 |   7.2367e-02 |       17 |   7.9325e+00
      94 |   2.8648e-07 |   3.3003e-07 |   8.0394e-01 |   1.0928e-01 |       16 |   7.4665e+00
      95 |   2.3034e-07 |   2.6536e-07 |   8.0404e-01 |   6.7303e-02 |       15 |   7.0007e+00
      96 |   1.8523e-07 |   2.1338e-07 |   8.0413e-01 |   6.4192e-02 |       14 |   6.5353e+00
      97 |   1.4896e-07 |   1.7161e-07 |   8.0421e-01 |   6.5703e-02 |       13 |   6.0700e+00
      98 |   1.1981e-07 |   1.3802e-07 |   8.0428e-01 |   9.6289e-02 |       12 |   5.6049e+00
      99 |   9.6365e-08 |   1.1101e-07 |   8.0434e-01 |   5.9033e-02 |       11 |   5.1400e+00
     100 |   7.7516e-08 |   8.9299e-08 |   8.0439e-01 |   6.4781e-02 |       10 |   4.6752e+00
     101 |   6.2357e-08 |   7.1836e-08 |   8.0444e-01 |   6.1663e-02 |        9 |   4.2105e+00
     102 |   5.0165e-08 |   5.7791e-08 |   8.0448e-01 |   6.6980e-02 |        8 |   3.7459e+00
     103 |   4.0358e-08 |   4.6494e-08 |   8.0452e-01 |   9.8091e-02 |        7 |   3.2815e+00
     104 |   3.2470e-08 |   3.7406e-08 |   8.0455e-01 |   6.7394e-02 |        6 |   2.8171e+00
     105 |   2.6125e-08 |   3.0096e-08 |   8.0457e-01 |   6.9847e-02 |        5 |   2.3527e+00
     106 |   2.1020e-08 |   2.4215e-08 |   8.0460e-01 |   7.3793e-02 |        4 |   1.8885e+00
     107 |   1.6913e-08 |   1.9484e-08 |   8.0461e-01 |   9.4014e-02 |        3 |   1.4243e+00
     108 |   1.3609e-08 |   1.5677e-08 |   8.0463e-01 |   5.8040e-02 |        2 |   9.6012e-01
     109 |   1.0950e-08 |   1.2615e-08 |   8.0465e-01 |   6.2136e-02 |        1 |   4.9600e-01
     110 |   8.8111e-09 |   1.0151e-08 |   8.0466e-01 |   6.4131e-02 |        0 |   3.1912e-02
------------------------------------------------------------------
Results of Time Iteration Algorithm
 * Complementarities: false
 * Discretized Process type: Dolo.DiscretizedIIDProcess
 * Decision Rule type: Dolo.CubicDR{Dolo.EmptyGrid{1}, Dolo.UCGrid{1}, 1, 1}
 * Number of iterations: 110
 * ϵₙ=|F(xₙ,xₙ)| < 1.0e-08: true
 * |x - x'| < 1.0e-08: false
 

```

$\lambda_n=\frac{\eta_n}{\eta_{n-1}}$ with $\eta_n = |{x_n-x_{n-1}}|$, so when lambda is equal or greater than one, in the current iteration we are further away from finding a fixed-point for the policy function $x_n$ than in the previous iteration. In this case, $\lambda_n$ is not a good indicator measure of convergence of policy function, as if we apply recursively $\eta_n = \lambda_n * \eta_{n-1}$, the sequence $(\eta_n)$ would be diverging. As a result, we cannot estimate the number of iterations before convergence of eta towards its tolerance level and thus we cannot estimate the remaining time. We change $p=NaN$ for lamba greater or equal to one.

We also change names of columns in the output table such that the row-name fits in one row.
 

 ```r
 
function time_iteration(model;
    dr0=Dolo.ConstantDecisionRule(model),
    discretization=Dict(),
    interpolation=:cubic,
    verbose=true,
    details=true,
    ignore_constraints=false,
    trace = false,
    tol_η = 1e-8,
    tol_ε = 1e-8,
    λbar = 8.0466e-01, 
    maxit=500
)
    
    F = Euler(model; discretization=discretization, interpolation=interpolation, dr0=dr0,  ignore_constraints=ignore_constraints)

    complementarities = false


    if interpolation != :cubic
        error("Interpolation option ($interpolation) is currently not recognized.")
    end

    ti_trace = trace ? IterationTrace([F.x0]) : nothing

    z0 = deepcopy(F.x0)
    z1 = deepcopy(z0)

    vector_time = rand(maxit) 

    local err_ε, err_η, z0, z1, p, it 

    log = Iterationlog(
        it = ("n", Int),
        err =  ("εₙ=|F(xₙ,xₙ)|", Float64),
        sa =  ("ηₙ=|xₙ-xₙ₋₁|", Float64),
        lam = ("λₙ=ηₙ/ηₙ₋₁", Float64),
        elapsed = ("Time (s)", Float64),
        nb_it_before_convergence_of_x = ("N-n", Int), 
        remaining_time = ("Remain (s)", Float64), 
    )

    initialize(log, verbose=verbose; message="Time Iteration")

    err_η_0 = NaN

    it = 0

    while it<=maxit

        it += 1

        t1 = time_ns()

        r0 = F(z0, z0; set_future=true)

        err_ε=norm(r0)

        fun = u->F(u, z0; set_future=false)
        dfun = u->df_A(F,u, z0; set_future=false)

        sol = newton2(
            fun, dfun,
            z0, verbose=false
        )

        z1 = sol
        δ = z1 - z0

        #trace && push!(ti_trace.trace, z1) #

        err_η = norm(δ)
        gain = err_η / err_η_0
        err_η_0 = err_η

        # z0.data[:] .= z1.data
        z0 = z1
        
        p = NaN

        if gain>=1
            p = NaN
        elseif gain<1
            p = Base.log(tol_η/err_η) / Base.log(λbar)
        end

        elapsed = time_ns() - t1

        elapsed /= 1000000000

        vector_time[it] = elapsed 

        local avg_time
        
        avg_time = mean(vector_time[max(1,end-4):end])

        time_left = NaN
        if gain>=1
            time_left = NaN
            p = NaN
        elseif gain<1
            time_left = avg_time*p
            p = round(Int, p) 
        end

     append!(log; 
        verbose=verbose, 
        it=it, 
        err=err_ε, 
        sa=err_η, 
        lam=gain, 
        elapsed=elapsed,
        nb_it_before_convergence_of_x=p, 
        remaining_time=time_left,
     )

     if err_ε<tol_ε     
         break
     end

     if err_η<tol_η       
         break 
     end

end

 finalize(log, verbose=verbose)

 res = TimeIterationResult(
     F.dr.dr, 
     it, 
     complementarities, 
     F.dprocess, 
     err_ε, 
     err_η, 
     tol_ε, 
     tol_η, 
     λbar, 
     p, 
     log, 
     ti_trace
     )

 return res

end
```