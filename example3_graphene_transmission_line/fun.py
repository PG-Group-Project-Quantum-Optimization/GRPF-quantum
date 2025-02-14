import numpy as np

def fun(z):
    """
    Compute the function value based on the given argument.

    Args:
    z (complex): Function argument.

    Returns:
    complex: Function value.
    """
    f = 1e12
    c = 299792458
    uo = 4 * np.pi * 1e-7
    eo = 1 / (uo * c * c)

    e = 1.602176565e-19
    kB = 1.3806488e-23
    hk = 1.05457168e-34
    vFe = 1e6
    muc = 0.05 * e
    t = 0.135e-12
    T = 300
    er1 = 1
    er2 = 11.9

    w = 2 * np.pi * f
    k0 = w / c
    kro = -1j * z * k0

    Slo = -1j * e * e * kB * T * np.log(2 + 2 * np.cosh(muc / kB / T)) / (np.pi * hk * hk * (w - 1j / t))

    a = -3 * vFe * vFe * Slo / (4 * (w - 1j / t)**2)
    b = a / 3

    Y1TM = w * er1 * eo / np.sqrt(er1 * k0**2 - kro**2)
    Y2TM = w * er2 * eo / np.sqrt(er2 * k0**2 - kro**2)
    YSTM = Slo + 1 * a * kro**2 + 1 * b * kro**2

    w = (Y1TM + Y2TM + YSTM) * (-Y1TM + Y2TM + YSTM) * (Y1TM - Y2TM + YSTM) * (-Y1TM - Y2TM + YSTM)

    return w