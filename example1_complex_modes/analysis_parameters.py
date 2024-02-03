import numpy as np
from disk_dom import disk_dom

# Circular disk domain definition
R = 1  # Domain radius
r = 0.15  # Initial mesh step

NewNodesCoord = disk_dom(R, r)  # Initial mesh generation

Tol = 1e-9  # Accuracy (candidate region size)

visual = 2  # Mesh visualization:  0 - turned off,   1  - only last iteration,   2  - all iterations

ItMax = 100  # Max number of iterations

NodesMax = 500000  # Max number of nodes

SkinnyTriangle = 3  # Skinny triangle definition
