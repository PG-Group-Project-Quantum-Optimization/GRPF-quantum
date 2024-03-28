def fun(z):
    """
    Compute the function value based on the given argument.

    Args:
    z (complex): Function argument.

    Returns:
    complex: Function value.

    Example:
    >>> result = fun(2 + 3j)
    """
    w = (z - 1) * (z - 1j) ** 2 * (z + 1) ** 3 / (z + 1j)
    return w