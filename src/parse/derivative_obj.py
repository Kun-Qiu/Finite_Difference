"""
Derivative object to represent the derivatives in any arbitary
linear partial differential equations.
"""


class Derivative:
    def __init__(self, deriv_tuple):
        self.coeff_sign, self.coeff_mag = deriv_tuple[0]
        self.diff_order = deriv_tuple[1]
        self.depen_var = deriv_tuple[2]
        self.indep_var = deriv_tuple[3]
        self.indep_var_order = deriv_tuple[4]
