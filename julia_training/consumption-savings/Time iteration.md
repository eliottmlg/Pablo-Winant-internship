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

Still the same problem, variable undefined. An explanation is that, in Julia, a variable x requires a new binding each time a loop is executed.


