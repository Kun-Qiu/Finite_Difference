import re
from src.parse.derivative_obj import Derivative

def extract_order(exp_str):
    """
    Extract the order of the derivative

    :param exp_str  :   Expression string
    :return         :   Order of derivative
    """
    if exp_str is None or exp_str.strip() == '':
        return 1
    return int(exp_str.strip('^'))


def split_sign_and_coefficient(term):
    """
    Split the term containing both sign and coefficient into separate parts.

    Args:
    term (tuple): A tuple where the first element contains the sign and coefficient.

    Returns:
    tuple: A new tuple with the sign and coefficient separated.
    """
    # Extract the part with the sign and coefficient
    sign_and_coeff, *rest = term

    # Match the sign and the coefficient separately
    match = re.findall(r'([+-]?)\s*([\w\d./*]*)', sign_and_coeff.strip())
    if not match or len(match[0]) < 2:
        raise ValueError("Invalid format for sign and coefficient.")

    # Extract the sign and coefficient
    sign = match[0][0] if match[0][0] else '+'
    coeff = match[0][1].rstrip('*') if match[0][1] else '1'
    return ([sign, coeff], *rest)


class Eqn_Parser:
    def __init__(self, linear_eqn):
        self.eqn = linear_eqn
        self.dt = None
        self.eqn_var = set()
        self.time_derivs = ([], [])
        self.x_derivs = ([], [])
        self.y_derivs = ([], [])
        self.z_derivs = ([], [])
        self.t = ([], [])
        self.x = ([], [])
        self.y = ([], [])
        self.z = ([], [])
        self.const = ([], [])


    def is_time_dependent(self):
        """
        Return: Boolean --> whether the pde is time dependent
        """

        return not(self.time_derivs[0] == [] or self.time_derivs[1] == [])


    def parse_derivative(self, str):
        """
        Parse the derivatives in a function and retain in the fashion of a tuple:
        tuple = (coefficient, exp_numerator, dependent_var, var, exp_denominator)

        :param str  :   Input string
        :Return     :   Parsed equation
        """
        pde_str = str.strip()

        # Pattern for any derivative: d(^n)?(var)/d(var)(^n)?
        deriv_pattern = r'([+-]?\s*[\w\d*/.]+)?\s*d(\^\d+)?([A-Za-z]+)\s*/\s*d([A-Za-z])(\^\d+)?'
        all_derivs = re.findall(deriv_pattern, pde_str)
        for i in range(len(all_derivs)):
            revised_deriv = (
            all_derivs[i][0], 
            extract_order(all_derivs[i][1]),  
            all_derivs[i][2],  
            all_derivs[i][3],  
            extract_order(all_derivs[i][4]) 
            )
            all_derivs[i] = revised_deriv
        
        # Dependent and Independent Variable in equation
        for dphi in all_derivs:
            self.eqn_var.update(dphi[2:4])
    
        non_derivative_str = re.sub(deriv_pattern, '', pde_str)
        non_derivative = re.findall(r'[+-]?\s*[\w\d*/.]+', non_derivative_str)

        # Process each non-derivative term
        non_derivs = []
        for term in non_derivative:
            # Match sign, coefficient, and variable
            match = re.match(r'([+-]?)\s*(\d*)([a-zA-Z]*)', term.strip())
            if not match:
                raise ValueError(f"Invalid term format: {term}")

            sign = match.group(1) if match.group(1) else '+'
            coeff = match.group(2) if match.group(2) else ''
            var = match.group(3) if match.group(3) else ''

            # Check if the variable is in self.eqn_var
            if var in self.eqn_var:
                non_derivs.append((sign, int(coeff) if coeff else 1, var))
            else:
                non_derivs.append((sign, f"{coeff}{var}", ''))

        return all_derivs, non_derivs


    def parse_pde(self):
        """
        Execution of the parsing function

        :param      :   None
        :return     :   None
        """

        pde_str = self.eqn.strip()
        
        lhs, rhs = pde_str.split('=', 1)
        lhs, rhs = lhs.strip(), rhs.strip()
    
        eqn_derivs_lhs, eqn_const_lhs = self.parse_derivative(lhs)
        eqn_derivs_rhs, eqn_const_rhs = self.parse_derivative(rhs)

        for i in range(len(eqn_derivs_lhs)):
            eqn_derivs_lhs[i] = split_sign_and_coefficient(eqn_derivs_lhs[i])
            if eqn_derivs_lhs[i][3] == 't':
                self.time_derivs[0].append(Derivative(eqn_derivs_lhs[i]))
            elif eqn_derivs_lhs[i][3] == 'x':
                self.x_derivs[0].append(Derivative(eqn_derivs_lhs[i]))
            elif eqn_derivs_lhs[i][3] == 'y':
                self.y_derivs[0].append(Derivative(eqn_derivs_lhs[i]))
            elif eqn_derivs_lhs[i][3] == 'z':
                self.z_derivs[0].append(Derivative(eqn_derivs_lhs[i]))

        for i in range(len(eqn_const_lhs)):
            if eqn_const_lhs[i][2] == 't':
                self.t[0].append(eqn_const_lhs[i])
            elif eqn_const_lhs[i][2] == 'x':
                self.x[0].append(eqn_const_lhs[i])
            elif eqn_const_lhs[i][2] == 'y':
                self.y[0].append(eqn_const_lhs[i])
            elif eqn_const_lhs[i][2] == 'z':
                self.z[0].append(eqn_const_lhs[i])
            else:
                self.const[0].append(eqn_const_lhs[i])
        
        for i in range(len(eqn_derivs_rhs)):
            eqn_derivs_rhs[i] = split_sign_and_coefficient(eqn_derivs_rhs[i])
            if eqn_derivs_rhs[i][3] == 't':
                self.time_derivs[1].append(Derivative(eqn_derivs_rhs[i]))
            elif eqn_derivs_rhs[i][3] == 'x':
                self.x_derivs[1].append(Derivative(eqn_derivs_rhs[i]))
            elif eqn_derivs_rhs[i][3] == 'y':
                self.y_derivs[1].append(Derivative(eqn_derivs_rhs[i]))
            elif eqn_derivs_rhs[i][3] == 'z':
                self.z_derivs[1].append(Derivative(eqn_derivs_rhs[i]))

        for i in range(len(eqn_const_rhs)):
            if eqn_const_rhs[i][2] == 't':
                self.t[1].append(eqn_const_rhs[i])
            elif eqn_const_rhs[i][2] == 'x':
                self.x[1].append(eqn_const_rhs[i])
            elif eqn_const_rhs[i][2] == 'y':
                self.y[1].append(eqn_const_rhs[i])
            elif eqn_const_rhs[i][2] == 'z':
                self.z[1].append(eqn_const_rhs[i])
            else:
                self.const[1].append(eqn_const_rhs[i])


    def print(self, eqn_dir='lhs'):
        """
        Print all the variables stored in the parsing object

        :param eqn_dir  :   Side of equation
        :return         :   Print all variable to console
        """    

        if eqn_dir == "lhs":
            dir = 0
        elif eqn_dir == "rhs":
            dir = 1
        else:
            raise ValueError("Invalid eqn_dir. Use 'lhs' or 'rhs'.")

        # Titles and variables to print
        variables = {
            "Time Derivatives": self.time_derivs,
            "X Derivatives": self.x_derivs,
            "Y Derivatives": self.y_derivs,
            "Z Derivatives": self.z_derivs,
            "Time Variables (t)": self.t,
            "X Variables (x)": self.x,
            "Y Variables (y)": self.y,
            "Z Variables (z)": self.z,
            "Constants": self.const,
        }

        for title, value in variables.items():
            print(f"{title} ({'LHS' if dir == 0 else 'RHS'}): {value[dir]}")



# ------------------ Examples ------------------

# Time-dependent PDE:
# pde1 = "d^2u/dt^2 + c*du/dx + 4u = du/dt + 4k"
# eqn_parser = Eqn_Parser(pde1)
# eqn_parser.parse_pde()
# eqn_parser.print('lhs')
# eqn_parser.print('rhs')
