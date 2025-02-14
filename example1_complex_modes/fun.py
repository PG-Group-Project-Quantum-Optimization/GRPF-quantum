import numpy as np
from scipy.special import jv, yv

def fun(z):
    """
    Compute the function value based on the given argument.

    Args:
    z (complex): Function argument.

    Returns:
    complex: Function value.

    """
    z = z * 10

    f = 5e9
    c = 3e8
    mu0 = 4 * np.pi * 1e-7
    eps0 = 1e-9 / (36 * np.pi)
    a = 6.35e-3
    b = 10.0e-3
    eps_r1 = 10
    eps_r2 = 1
    m = 1

    omega = 2 * np.pi * f
    k0 = omega / c
    alpha = np.real(z) * k0
    beta = np.imag(z) * k0
    gamma = alpha + 1j * beta
    eps1 = eps0 * eps_r1
    eps2 = eps0 * eps_r2
    mu1 = mu0
    mu2 = mu0
    kappa1 = np.sqrt(gamma**2 + k0**2 * eps_r1)
    kappa2 = np.sqrt(gamma**2 + k0**2 * eps_r2)
    eta1 = np.sqrt(mu1 / eps1)
    eta2 = np.sqrt(mu2 / eps2)

    Jm_a1 = jv(m, kappa1 * a)
    Jm_a2 = jv(m, kappa2 * a)
    Ym_a2 = yv(m, kappa2 * a)
    Jm_b2 = jv(m, kappa2 * b)
    Ym_b2 = yv(m, kappa2 * b)
    DJm_a1 = (jv(m - 1, kappa1 * a) - jv(m + 1, kappa1 * a)) / 2
    DJm_a2 = (jv(m - 1, kappa2 * a) - jv(m + 1, kappa2 * a)) / 2
    DJm_b2 = (jv(m - 1, kappa2 * b) - jv(m + 1, kappa2 * b)) / 2
    DYm_a2 = (yv(m - 1, kappa2 * a) - yv(m + 1, kappa2 * a)) / 2
    DYm_b2 = (yv(m - 1, kappa2 * b) - yv(m + 1, kappa2 * b)) / 2

    W = np.array([
        [Jm_a1, 0, -Jm_a2, -Ym_a2, 0, 0],
        [0, Jm_a1 / eta1, 0, 0, -Jm_a2 / eta2, -Ym_a2 / eta2],
        [gamma * m * Jm_a1 / (a * kappa1**2), -omega * mu1 * DJm_a1 / (kappa1 * eta1),
         -gamma * m * Jm_a2 / (a * kappa2**2), -gamma * m * Ym_a2 / (a * kappa2**2),
         omega * mu2 * DJm_a2 / (kappa2 * eta2), omega * mu2 * DYm_a2 / (kappa2 * eta2)],
        [-omega * eps1 * DJm_a1 / kappa1, -m * gamma * Jm_a1 / (a * kappa1**2 * eta1),
         omega * eps2 * DJm_a2 / kappa2, omega * eps2 * DYm_a2 / kappa2,
         m * gamma * Jm_a2 / (a * kappa2**2 * eta2), m * gamma * Ym_a2 / (a * kappa2**2 * eta2)],
        [0, 0, Jm_b2, Ym_b2, 0, 0],
        [0, 0, gamma * m * Jm_b2 / (b * kappa2**2), gamma * m * Ym_b2 / (b * kappa2**2),
         -omega * mu2 * DJm_b2 / (kappa2 * eta2), -omega * mu2 * DYm_b2 / (kappa2 * eta2)]
    ])

    w = np.linalg.det(W)
    return w