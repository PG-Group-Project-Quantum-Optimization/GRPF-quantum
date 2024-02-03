import numpy as np
from rect_dom import rect_dom

# Bounding coordinates of analysed area.
xb = -2
xe = 2
yb = -2
ye = 2


# Initial distance between nodes
#r = 0.2

# DEBUG
r = 1.5
# DEBUG

NewNodesCoord = rect_dom(xb, xe, yb, ye, r)

Tol = 1e-9
visual = 2
ItMax = 10
NodesMax = 500000
SkinnyTriangle = 3
