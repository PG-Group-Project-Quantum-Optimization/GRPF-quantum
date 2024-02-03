import numpy as np
from rect_dom import rect_dom

# Bounding coordinates of the analyzed area.
xb = -100
xe = 400
yb = -100
ye = 400

# Initial distance between nodes
r = 18

NewNodesCoord = rect_dom(xb, xe, yb, ye, r)

Tol = 1e-9
visual = 2
ItMax = 100
NodesMax = 500000
SkinnyTriangle = 3
