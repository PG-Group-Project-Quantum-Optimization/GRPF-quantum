import numpy as np

def fun(z):
    """
    Compute the function value based on the given argument.

    Args:
    z (complex): Function argument.

    Returns:
    complex: Function value.
    """
    w = ((z - (4.0+3.0j))**2)*(z + 3.0)*((z + 1.0j)**3)*((z - 2.0)**2)*(z + (-2.0-1.0j))
    return w