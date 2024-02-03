def vinq(f):
    """
    Converts the function values to quadrants.

    Args:
        f (float): Function value.

    Returns:
        int: Quadrant.
    """

    if f.real > 0 and f.imag >= 0:
        w = 1
    elif f.real <= 0 and f.imag > 0:
        w = 2
    elif f.real < 0 and f.imag <= 0:
        w = 3
    elif f.real >= 0 and f.imag < 0:
        w = 4
    else:
        w = None

    return w
