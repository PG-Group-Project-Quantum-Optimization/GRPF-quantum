import numpy as np
from rect_dom import rect_dom  # Assuming rect_dom is the Python equivalent of the MATLAB function

# Bounding coordinates of the analyzed area.
xb = -2
xe = 2
yb = -2
ye = 2

# Initial distance between nodes
r = 0.1

NewNodesCoord = rect_dom(xb, xe, yb, ye, r)  # Initial mesh generation

Tol = 1e-9  # Accuracy
visual = 2  # Mesh visualization: 0 - turned off, 1 - only last iteration, 2 - all iterations
ItMax = 100  # Max number of iterations
NodesMax = 500000  # Max number of nodes
SkinnyTriangle = 3  # Skinny triangle definition
