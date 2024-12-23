import matplotlib.pyplot as plt
import numpy as np
# from src.solver.inheritance.heat_eqn_1D import HeatEqn1D

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

# heat_eqn_1d = HeatEqn1D(k, rho, c_p, N, dt, time, length_1d, T_initial_1d)
# heat_eqn_1d.solve(method='semi')

# length_2d = [2.0, 1.0]  # Length of the domain
# T_initial_2d = np.zeros((N, N))
# T_initial_2d[0, :] = T_initial_2d[-1, :] = T_initial_2d[:, 0] = T_initial_2d[:, -1] = init_temp

# heat_eqn_2d = HeatEqn2D(k, rho, c_p, N, dt, time, length_2d, T_initial_2d)
# heat_eqn_2d.solve()

from src.parse.derivative_obj import Derivative
from src.discretize.discretize_scheme.forward_difference import forward_diff

# pde1 = "du/dt + adu/dx = 0"

time_der = Derivative((['+', '1'], '1', 'u', 't', '1'))
space_der = Derivative((['+', 'a'], '1', 'u', 'x', '1'))
rhs_expr = 0 

lhs_expr = forward_diff(time_der)
print(rhs_expr - forward_diff(space_der))

