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







