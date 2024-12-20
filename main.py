import matplotlib.pyplot as plt
import numpy as np
from src.inheritance.heat_eqn_1D import HeatEqn1D

"""
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
"""

# Example Usage
k = 237
rho = 2710
c_p = 900
N = 100
dt = 0.1
time = 2000
init_temp = 100
length_1d = [1.0]
T_initial_1d = np.zeros(N)
T_initial_1d[0] = T_initial_1d[-1] = init_temp

heat_eqn_1d = HeatEqn1D(k, rho, c_p, N, dt, time, length_1d, T_initial_1d)
heat_eqn_1d.solve(method='semi')

# length_2d = [2.0, 1.0]  # Length of the domain
# T_initial_2d = np.zeros((N, N))
# T_initial_2d[0, :] = T_initial_2d[-1, :] = T_initial_2d[:, 0] = T_initial_2d[:, -1] = init_temp

# heat_eqn_2d = HeatEqn2D(k, rho, c_p, N, dt, time, length_2d, T_initial_2d)
# heat_eqn_2d.solve()
