from src.solver.eqn_base import HeatEqnBase
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


class HeatEqn2D(HeatEqnBase):
    def __init__(self, k, rho, c_p, N, dt, t, length, T_i):
        super().__init__(k, rho, c_p, N, dt, t, length, T_i)
        self.SPACE_STEP_X_1 = self.LENGTH[0] / (self.NUM_PT - 1)
        self.SPACE_STEP_X_2 = self.LENGTH[1] / (self.NUM_PT - 1)
        self.colorbar = None

    def construct_grid(self):
        x = np.linspace(0, self.LENGTH[0], self.NUM_PT)
        y = np.linspace(0, self.LENGTH[1], self.NUM_PT)
        self.GRID = np.meshgrid(x, y, indexing='ij')
        return self.GRID

    def build_matrix(self):
        coeffs = (self.ALPHA * self.TIME_STEP) * \
                 np.array([1 / np.square(self.SPACE_STEP_X_1),
                           1 / np.square(self.SPACE_STEP_X_2)])

        grid_size = np.power(self.NUM_PT, self.DIM)
        offsets = [0, -1, 1, -self.NUM_PT, self.NUM_PT]

        def general_matrix(main_diagonal, off_diagonal_x, off_diagonal_y, NUM_PT):
            diag = np.ones(grid_size)
            off_diag_1 = np.zeros(grid_size)
            off_diag_2 = np.zeros(grid_size)

            non_boundary_indices = np.where(
                (np.arange(grid_size) % NUM_PT != 0) &  # Not on left boundary
                ((np.arange(grid_size) + 1) % NUM_PT != 0) &  # Not on right boundary
                (np.arange(grid_size) >= NUM_PT) &  # Not in top boundary row
                (np.arange(grid_size) < grid_size - NUM_PT)  # Not in bottom boundary row
            )

            diag[non_boundary_indices] *= main_diagonal
            off_diag_1[non_boundary_indices] = off_diagonal_x
            off_diag_2[non_boundary_indices] = off_diagonal_y

            return [diag, off_diag_1[1:], off_diag_1[:-1],
                    off_diag_2[self.NUM_PT:], off_diag_2[:(-self.NUM_PT)]]

        self.A = diags(general_matrix((1 + 2 * (coeffs[0] + coeffs[1])),
                                      (-coeffs[0]), -coeffs[1], self.NUM_PT),
                       offsets=offsets, format='csc')
        self.Ac = diags(general_matrix((1 - 2 * (coeffs[0] + coeffs[1])),
                                       (coeffs[0]), coeffs[1], self.NUM_PT),
                        offsets=offsets, format='csc')

    def solve(self):
        plt.ion()
        self.construct_grid()  # Ensure grid is initialized
        self.build_matrix()

        mat_shape = np.shape(self.b)
        while self.time < self.TIME:
            rhs = self.Ac.dot(self.b.flatten())
            T_new = spsolve(self.A, rhs)
            self.b = T_new.reshape(mat_shape)

            self.update_plot(self.time)
            plt.pause(0.1)

            self.time += self.TIME_STEP

    def update_plot(self, time):
        """Updates the plot with the current temperature profile and time."""
        self.ax.clear()
        contour = self.ax.contourf(self.GRID[0], self.GRID[1], self.b, cmap='hot')

        # Remove previous colorbar if it exists
        if self.colorbar is not None:
            self.colorbar.remove()

        # Add new colorbar and store the reference
        self.colorbar = plt.colorbar(contour, ax=self.ax, label="Temperature")

        self.ax.set_title(f"Temperature at t={time:.2f}")
        self.ax.set_xlabel("X Position")
        self.ax.set_ylabel("Y Position")
        plt.draw()