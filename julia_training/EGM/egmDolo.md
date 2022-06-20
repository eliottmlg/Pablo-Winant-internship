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

The EGM is a solution method for a certain class of models. 

Consider a model of the following form:

$z_{t} = E_t[h(m_{t+1},s_{t+1},x_{t+1})]$, the expectation equation

$s_t = g_t(m_{t-1},a_{t-1},m_t)$, the state half-transition equation

$x_t = \tau(m_t,a_t,z_t)$, the direct-response equation

$s_t = a\tau(m_t,a_t,x_t)$, the reverse-state equation

where $m_t$ is the set of exogenous shocks at t, $x_t$ is the set of controls at t, $s_t$ is the set of states.

The EGM consists in approximating policy functions (optimal controls) as a function of exogenous and endogenous controls, $x_t = \phi(m_t,s_t)$, alike the Time iteration algorithm. The difference lies the fact that optimal decisions are computed only at gridpoints that satisfy the Euler equation and these computations do not involve any root-finding problem. At step $n$, the current guess for the control, $x(s_t)=\phi^n(m_t,s_t)$, serves as the control being used in the next period.

Here is an outline of the algorithm:

1. Set a fixed grid over the continuous post-state variable $a_t$, say $A$
2. Set with an initial guess of the policy function $\phi^0$, a mapping from the grid of the state to the set of decision variable values 
3. Compute at the current guess $\phi^n$, for each point on the fixed grid $A$, the RHS of the Euler equation, $z_t = E_t[h(m_{t+1},s_{t+1},x_{t+1})]$ 
4. Compute the updated decision at each gridpoint of $A$ using the direct-response equation
5. Construct the new grid of the state $s_t$ using the reverse-state equation
6. Match each endogenous grid points of $s_t$ to decision values $x_t$ and interpolate the policy function $x_n(s_t)$
7. Check for convergence, i.e. if $\eta_n = |x_n(s_t) - x_{n-1}(s_t)|$ is below tolerance level, if so stop, if not start over from 3. 