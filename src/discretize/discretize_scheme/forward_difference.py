from sympy import symbols

def forward_diff(deriv):
    """
    Calculate the forward difference approximation of a derivative.
    
    :param deriv    : A derivative object
    :return         : A sympy expression for the forward difference.
    """
    
    coeff, var_delta, var_0, step = symbols(
        f"{deriv.coeff_mag}, {deriv.depen_var}_(i+1), {deriv.depen_var}_i, d{deriv.indep_var}")
    
    forward_diff_expr = (coeff / step) * (var_delta - var_0) 

    return forward_diff_expr


