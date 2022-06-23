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

$s_t = g(m_{t-1},a_{t-1},m_t)$, the state half-transition equation **# why is a indexed by t?**

$x_t = \tau(m_t,a_t,z_t)$, the direct-response equation

$s_t = a\tau(m_t,a_t,x_t)$, the reverse-state equation

where $m_t$ is a vector of length $K$ of exogenous shocks at t, $x_t$ is the control at t, $s_t$ is the state at t, $a$ is the continuous post-state.

The EGM consists in approximating the policy function as a function of shocks and the state, $\phi(m_t,s_t)$, alike the Time iteration algorithm. The difference lies in the fact that optimal decisions are computed only at gridpoints that satisfy the Euler equation and these computations do not involve any root-finding problem. Decisions can be written in two equivalent ways, in terms of the state or the post-state $x^{i,j}=\phi^n(m_j,s^i) = \psi^n(m_j,a^i)$. Here, $j=1,2,...,J$ indicates the current state of the shock, $k=1,2,...,K$ indicates next period's state of the shock, and $i=1,2,...,I$ indicates the point on the fixed grid $A$. At step $n$, the current guess for the control, $\phi_n(m^j,s^i)$, serves as the control being used in the next period.
Here is an outline of the algorithm:

1. Discretize the post-state variable $a$ on a fixed grid $A$, the state $s_t^0$ on an endogenous grid $S^0$, and the shocks $m$ on exogenous grid $M$
2. Start with an initial guess of the policy function $\phi^0(m_t,s_t)$, a mapping from the grid of the state and the exogenous variable to the set of decision variable values 
3. At the current guess $\phi_n(m,s)$, execute the following for $j=1$ through $j=J$: $\forall a^i \in A$, compute $z(m_j,a^i) = E_t[h(m_{t+1},s_{t+1},x_{t+1})] = E_{m_k|m_j}[h(m_k,g(m_j,a^i,m_k),\phi^{n}(m_k,g(m_j,a,m_k)))]$ for $k=1,...,K$. This corresponds to the RHS of the Euler equation
4. Compute for all $(a^i, m_j)$ the updated decisions $x_{n+1}^{i,j}$ using the direct-response equation: $x_{n+1}(m_j, a^i) = \tau(m_j,a^i,z(m_j,a^i))$
5. Construct the new endogenous grid $S^{n+1}$ by computing for all $(m_j,a^i)$ the state $s_t^{n+1}$ using the reverse-state equation: $s_t^{n+1} = a\tau(m_t,a,x_{n+1}^{i,j})$
6. Check for convergence, i.e. if $\eta_n = |x_{n+1}(s_t) - x_{n}(s_t)|$ is below tolerance level. If so, stop the loop and match each endogenous gridpoint of $s_t$ with an element of the decision vector $x_t$, finally use interpolation to construct the policy function $x_n(s_t)$. If not, start over from 3. 