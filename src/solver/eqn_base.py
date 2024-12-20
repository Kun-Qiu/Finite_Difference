import matplotlib.pyplot as plt
import numpy as np

# Virtual class for the construction of 1D to 3D solver for equations

class HeatEqnBase:
    def __init__(self, k, rho, c_p, N, dt, t, length, T_i):
        if isinstance(length, (float, int)):
            self.DIM = 1
            self.LENGTH = [length]
        elif hasattr(length, "__iter__"):
            self.DIM = len(length)
            self.LENGTH = length

        self.NUM_PT = N                 # Number of points
        self.TIME_STEP = dt             # Time step
        self.TIME = t                   # Ending Criteroia: time
        self.ALPHA = k / (rho * c_p)    # Thermal coefficient
        self.b = T_i                    # Initial temperature distribution
        self.A = None                   # Implicit -> lhs matrix
        self.Ac = None                  # Semi Implicit -> rhs matrix
        self.LIMIT_Y = np.max(T_i)      # Maximum Temperature
        self.GRID = None                # Grid structure
        self.time = 0                   # Current time
        
        self.fig, self.ax = plt.subplots()


    def construct_grid(self):
        raise NotImplementedError("This method should be implemented in child classes.")


    def build_matrix(self):
        raise NotImplementedError("This method should be implemented in child classes.")


    def solve(self):
        raise NotImplementedError("This method should be implemented in child classes.")


    def update_plot(self, time):
        raise NotImplementedError("This method should be implemented in child classes.")