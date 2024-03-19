# GRPF-quantum
A repository for code related to GRPF (Global Root and Pole Finding) algorithm

## How to use
### Setting up a complex function
To define a function that can be used in the code you need to create `fun()` function that takes in a complex number and outputs a complex number and a file **analysis_parameters.py** that contains information on the parameters for the analysis of `fun()` in the following format:
```python
import numpy as np
from rect_dom import rect_dom

# Bounding coordinates of the analyzed area
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
```

In the `grpf()` function import `ItMax`, `NodesMax`, `Tol`, `SkinnyTriangle`, `visual`, `NewNodesCoord` parameters and `fun()` function.

### Running the code
To run the code on your complex function simply run the **grpf.py** file.