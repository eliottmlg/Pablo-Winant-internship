(to open preview window: ctrl+shift+v)

# Implementing EGM in Dolo

The Dolo program reads YAML files written with the following structure:

```r
symbols:
exogenous: []
states: []
controls: []
parameters: []

equations:
transitions: 
-...
arbitrage: 
-...
half_transition: |
...
reverse_state: |
...
expectation: |
...
direct_response_egm: |
...

calibration:

domain:

exogenous: 

options:
```

## How is Dolo used to solve a model in Julia?
The modeller enters the equations describing the model and the values of parameters in the dedicated section, after which she uses an algorithm to solve the model. Algorithms include Value iteration (VI), Time iteration (TI), or Endogenous gridpoint method (EGM). To start the resolution, Dolo first needs to reads the model. It does so using a compiler-file, which extracts the equations and the parameter values from the YAML file. This information will then be called by the algorithm-file, which is the file describing the method of resolution chosen by the modeller (EGM, TI, VI).
The algorithm-file will run resolution routines and will print the results in the Julia REPL.

## EGM

The EGM is a solution method for a certain class of models. In particular, it handles models with one continuous state variable and one continous control variable.

Consider a model of the following form:

$z_{t} = E_t[h(m_{t+1},s_{t+1},x_{t+1})]$, the expectation equation

$s_t = g_t(m_{t-1},a_{t-1},m_t)$, the state half-transition equation **# why is a indexed by t?**

$x_t = \tau(m_t,a_t,z_t)$, the direct-response equation

$s_t = a\tau(m_t,a_t,x_t)$, the reverse-state equation

where $m_t$ is a vector of length $K$ of exogenous shocks at t, $x_t$ is the control at t, $s_t$ is the state at t, $a$ is the continuous post-state.

The EGM consists in approximating policy functions (optimal controls) as a function of exogenous and endogenous controls, $x_t = \phi(m_t,s_t)$, alike the Time iteration algorithm. The difference lies the fact that optimal decisions are computed only at gridpoints that satisfy the Euler equation and these computations do not involve any root-finding problem. At step $n$, the current guess for the control, $x(s_t)=\phi^n(m_t,s_t)$, serves as the control being used in the next period.

Here is an outline of the algorithm:

1. Discretize the post-state variable $a$ on a fixed grid $A$, the state $s_t^0$ on an endogenous grid $S^0$, and the shocks $m_t^1, m_t^2, ...,m_t^K$ on exogenous grids $M_1,M_2, ..., M_K$
2. Start with an initial guess of the policy function $\phi^0(m_t,s_t)$, a mapping from the grid of the state and the exogenous variable to the set of decision variable values 
3. At the current guess $\phi^n(m_t,s_t)$, execute the following: $\forall a \in A, m_t^k \in M_k, m_{t+1}^k \in M_k, k=1,...,K$, compute $z(m_t,a_t) = E_t[h(m_{t+1},s_{t+1},x_{t+1})] = E_{m_{t+1}|m_t}[h(m_{t+1},g_t(m_{t},a_{t},m_{t+1}),\phi^{n}(m_{t+1},g_t(m_{t},a_{t},m_{t+1})))]$, which corresponds to the RHS of the Euler equation
4. Compute for all $(a_t, m_t)$ the updated decisions $x_t^{n+1} = \phi^{n+1}(m_t,s_t)$ using the direct-response equation: $x^{n+1}(m_t, a_t) = \tau(m_t,a_t,z(m_t,a_t))$
5. Construct the new endogenous grid $S^{n+1}$ by computing for all $(a, m_t)$ the state $s_t^{n+1}$ using the reverse-state equation: $s_t^{n+1} = a\tau(m_t,a_t,x_t^{n+1})$
6. Check for convergence, i.e. if $\eta_n = |x^{n+1}(s_t) - x^{n}(s_t)|$ is below tolerance level. If so, stop the loop and match each endogenous gridpoint of $s_t$ with an element of the decision vector $x_t$, finally use interpolation to construct the policy function $x_n(s_t)$. If not, start over from 3. 