### Finite_Difference
Method of finite difference for solving linear ordinary differential equation / partial differential equations.


## Heat Equation
Heat equation subjected to the following

If using the explicit method check stability condition for explicit method 
(Courant–Friedrichs–Lewy condition)

The Crank Nicolson approximate solutions can still contain spurious oscillations 
if the ratio of (time step Δt * thermal diffusivity) over the square of space step 
Δx^2 is larger than 1/2

Dirichlet Boundary Conditions:
u[t, 0] = 0 = u[t, length]

Neumann Boundary Conditons:
dT/dt[t, 0] = 0 = dT/dt[t, length]

