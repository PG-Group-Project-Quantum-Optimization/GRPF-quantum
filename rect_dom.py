import numpy as np

def rect_dom(xb, xe, yb, ye, r):
    '''
    rect_dom: generates the initial mesh for rectangular domain z=x+jy x\in[xb,xe] , y\in[yb,ye]
    
    Args:
        xb     (float): real part range begin 
        xe     (float): real part range end 
        yb     (float): imag part range begin 
        ye     (float): imag part range end 
        r      (float): initial mesh step
    
    Returns:
        NewNodesCoord [(float, float)]: generated nodes coordinates
    '''
    X = xe - xb
    Y = ye - yb

    n = np.ceil(Y / r + 1)
    dy = Y / (n - 1)
    m = np.ceil(X / np.sqrt(r**2 - dy**2 / 4) + 1)
    dx = X / (m - 1)

    vx = np.linspace(xb, xe, int(m))
    vy = np.linspace(yb, ye, int(n))
    x, y = np.meshgrid(vx, vy)

    temp = np.ones((int(n), 1))
    temp[-1] = 0
    y = y + 0.5 * dy * np.kron((1 + (-1)**np.arange(1, m + 1)) / 2, temp)

    x = np.reshape(x, (int(m * n), 1))
    y = np.reshape(y, (int(m * n), 1))
    tx = ((np.arange(2, m + 1, 2) - 1) * dx + xb).reshape(-1, 1)
    ty = np.zeros_like(tx) + yb
    x = np.vstack([x, tx])
    y = np.vstack([y, ty])

    NewNodesCoord = np.concatenate([x, y], axis=1)

    return NewNodesCoord
