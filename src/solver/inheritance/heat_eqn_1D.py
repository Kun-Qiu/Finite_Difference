from src import eqn_base as eqn
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


class HeatEqn1D(eqn.HeatEqnBase):
    def __init__(self, k, rho, c_p, N, dt, t, length, T_i):
        super().__init__(k, rho, c_p, N, dt, t, length, T_i)
        self.SPACE_STEP_X_1 = self.LENGTH[0] / (self.NUM_PT - 1)

    def construct_grid(self):
        grid_range = np.linspace(0, self.LENGTH[0], self.NUM_PT)
        self.GRID = grid_range

        return grid_range

    def build_matrix(self, method='semi'):
        coeffs = (self.ALPHA * self.TIME_STEP) / np.square(self.SPACE_STEP_X_1)
        grid_size = np.power(self.NUM_PT, self.DIM)
        offsets = [0, -1, 1]

        def general_matrix(main_diagonal, off_diagonal):
            # Initialize main diagonal and set values based on main_diagonal parameter
            diag = np.ones(grid_size)
            diag[1:self.NUM_PT - 1] *= main_diagonal

            # Initialize off-diagonal values
            off_diag = np.ones(grid_size - 2) * off_diagonal

            return [diag, np.append(off_diag, 0), np.append(0, off_diag)]

        if method == 'im' or method == 'semi':
            self.A = diags(general_matrix((1 + coeffs), (-coeffs / 2)), offsets=offsets, format='csc')
            if method == 'semi':
                self.Ac = diags(general_matrix((1 - coeffs), (coeffs / 2)), offsets=offsets, format='csc')


    def solve(self, method='ex'):    
        """
        'ex' represents explicit solver and 'im' represents implicit and 'semi' represent Crank Nicolson
        """
        
        self.construct_grid()
        grid_space = self.GRID[2:] - self.GRID[1:-1]

        if (method == 'im') or (method == 'semi'):
            self.build_matrix(method=method) 
        
        plt.ion()

        while self.time < self.TIME:
            if method == 'ex':
                # Explicit Solver
                sol = self.b
                sol[1:-1] = self.ALPHA * self.TIME_STEP *  (sol[2:] - 2 * sol[1:-1] 
                                                            + sol[0:-2]) / (grid_space ** 2) + sol[1:-1]
                self.b = sol

            if method == 'im':
                # Implicit Solver
                rhs = self.b
                sol = spsolve(self.A, rhs)
                self.b = sol

            if method == 'semi':
                # Semi Implicit Solver
                rhs = self.Ac.dot(self.b)
                sol = spsolve(self.A, rhs)
                self.b = sol

            self.update_plot(self.time)
            plt.pause(0.1)

            self.time += self.TIME_STEP


    def update_plot(self, time):
        """Updates the plot with the current temperature profile and time."""
        self.ax.clear()
        self.ax.plot(self.GRID, self.b, label=f"t={time:.2f}")
        self.ax.set_xlabel("Position")
        self.ax.set_ylabel("Temperature")
        self.ax.set_ylim(0, self.LIMIT_Y)  # Set y-axis limits
        self.ax.legend()
        plt.draw()