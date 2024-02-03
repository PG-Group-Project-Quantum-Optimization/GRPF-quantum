import numpy as np

def fun(z):
    """
    Compute the function value based on the given argument.

    Args:
    z (complex): Function argument.

    Returns:
    complex: Function value.
    """
    ns = 0.065 - 4j
    n1 = 1.5835
    nc = 1.0
    d1 = 1.81e-6
    lambda0 = 0.6328e-6
    k0 = 2. * np.pi / lambda0
    k0d1 = k0 * d1
    kappa1 = np.sqrt(n1**2. - z**2.)
    gammas = np.sqrt(z**2. - ns**2.)
    gammac = np.sqrt(z**2. - nc**2.)
    m11 = np.cos(kappa1 * k0d1)
    m12 = 1j / kappa1 * np.sin(kappa1 * k0d1)
    m21 = 1j * kappa1 * np.sin(kappa1 * k0d1)
    m22 = np.cos(kappa1 * k0d1)

    w = np.linalg.det([[1.0, -m11 + 1j * gammac * m12], [1j * gammas, -m21 + 1j * gammac * m22]])

    return w
