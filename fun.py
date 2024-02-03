import numpy as np

def fun(z):
    # function definition

    f = 1.0e9
    er = 5 - 2j
    mr = 1 - 2j
    d = 1.0e-2
    c = 3e8
    w = 2 * np.pi * f
    k0 = w / c
    c = er**2 * (k0 * d)**2 * (er * mr - 1)
    w = er**2 * z**2 + z**2 * np.tan(z)**2 - c

    return w