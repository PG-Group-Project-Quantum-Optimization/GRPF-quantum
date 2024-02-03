import numpy as np

def find_next_node(nodes_coord, prev_node, ref_node, temp_nodes):
    """
    Finds the next node in the candidate region boundary process.
    The next node (after the reference one) is picked from the fixed set of nodes.

    Args:
        nodes_coord (list): Coordinates of nodes.
        prev_node (int): Previous node.
        ref_node (int): Reference (current) node.
        temp_nodes (set): Set of nodes.

    Returns:
        int: Index of the next node.
    """

    p = nodes_coord[prev_node, :]
    s = nodes_coord[ref_node, :]
    n = nodes_coord[temp_nodes, :]

    no_of_temp_nodes = n.shape[0]

    sp = np.ones((no_of_temp_nodes, 1)) * (p - s)
    sn = n - np.ones((no_of_temp_nodes, 1)) * s

    len_sp = np.sqrt(sp[:, 0]**2 + sp[:, 1]**2)
    len_sn = np.sqrt(sn[:, 0]**2 + sn[:, 1]**2)

    dot_prod = sp[:, 0] * sn[:, 0] + sp[:, 1] * sn[:, 1]

    phi = np.arccos(dot_prod / (len_sp * len_sn))

    temp = np.where(sp[:, 0] * sn[:, 1] - sp[:, 1] * sn[:, 0] < 0)
    phi[temp] = 2 * np.pi - phi[temp]

    index = np.argmin(phi)

    return index
