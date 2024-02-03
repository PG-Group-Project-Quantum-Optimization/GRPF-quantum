import numpy as np

def disk_dom(R, r):
    """
    Generates the initial mesh for a circular disk domain of radius R centered at (0,0).

    The mesh step size is given by 'r'. The generated nodes are within the region |z| < R.

    Args:
        R (float): Domain radius.
        r (float): Initial mesh step.

    Returns:
        list: Generated nodes' coordinates.
    """

    h = r * np.sqrt(3) / 2
    n = 1 + round(R / h)
    Rn = (np.arange(1, n + 1) * R / n).reshape(-1, 1)
    NewNodesCoord = np.array([[0, 0]])
    f0 = 0
    npoints = 6

    for i in range(n):
        f = f0 + np.linspace(0, 2 * np.pi, npoints + 1)[:-1]
        xyn = Rn[i] * np.column_stack((np.cos(f), np.sin(f)))
        NewNodesCoord = np.concatenate((NewNodesCoord, xyn))
        f0 = f0 + np.pi / (6 * n)
        npoints = npoints + 6

    return NewNodesCoord
