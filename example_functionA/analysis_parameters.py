import numpy as np
from rect_dom import rect_dom

# Bounding coordinates of analysed area.
xb = -8
xe = 8
yb = -8
ye = 8

r = 2

NewNodesCoord = rect_dom(xb, xe, yb, ye, r)

Tol = 1e-9
visual = 2
ItMax = 1000
NodesMax = 500000
SkinnyTriangle = 3
