import numpy as np
from rect_dom import rect_dom

# Bounding coordinates of the analyzed area.
xb = 1
xe = 2.5
yb = -1
ye = 1

# Initial distance between nodes
r = 0.5

NewNodesCoord = rect_dom(xb, xe, yb, ye, r)

Tol = 1e-9
visual = 2
ItMax = 1000
NodesMax = 500000
SkinnyTriangle = 3
