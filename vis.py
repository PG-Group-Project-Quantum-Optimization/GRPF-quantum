import numpy as np
import matplotlib.pyplot as plt


def vis(NodesCoord, Edges, Quadrants, PhasesDiff, it, Show=True, QuantumCandidate=None):
    '''
    Visualises the existing mesh, with colored edges representing quadrants and candidate edges.

    Args:
        NodesCoord [(float, float)]: Coordinates of the nodes of analysed mesh
        Edges [(int, int)]: Indexes of nodes connected with edges
        Quadrants [int]: Quadrants of value of the function for nodes in NodesCoord
        PhasesDiff [int]: Difference of phase between coordinates at the end of the edges described in Edges 
        QuantumCandidate: Edges candidate from quantum circuit.
        Show (bool): Whether to show the visualisation or wait for matplotlib.pyplot.show() after the function
    '''
    
    NoOfEdges = Edges.shape[0]
    EdgesColor = np.zeros((NoOfEdges, 1))

    if QuantumCandidate is not None:
        for i in QuantumCandidate:
            EdgesColor[
                np.where((Edges == np.sort(i)).all(1))
            ] = 5  # suspected edges
    else:
        EdgesColor[(PhasesDiff == 2) | np.isnan(PhasesDiff)] = 5  # suspected edges


    Temp = (PhasesDiff == 0)
    EdgesColor[Temp.flatten()] = Quadrants[Edges[Temp.flatten(), 0]]

    vNaN = np.full(NoOfEdges, np.nan)
    
    EdgesXCoord = np.column_stack((NodesCoord[Edges[:, 0], 0], NodesCoord[Edges[:, 1], 0], vNaN))
    EdgesYCoord = np.column_stack((NodesCoord[Edges[:, 0], 1], NodesCoord[Edges[:, 1], 1], vNaN))
   
    EdgesXCoord = EdgesXCoord.reshape((3 * NoOfEdges, 1))
    EdgesYCoord = EdgesYCoord.reshape((3 * NoOfEdges, 1))

    
    EdgesColor = np.kron(EdgesColor, np.array([[1], [1], [1]]))

    plt.figure(1)
    plt.clf()

    def subplot(number, color):
        plt.plot(EdgesXCoord[(EdgesColor == number) | (EdgesColor == np.nan)], EdgesYCoord[(EdgesColor == number) | (EdgesColor == np.nan)], color)
        
    subplot(0, 'k')
    subplot(1, 'r')
    subplot(2, 'y')
    subplot(3, 'g')
    subplot(4, 'b')
    subplot(5, 'm')
    plt.xlabel("Re(z)")
    plt.ylabel("Im(z)")
    plt.tight_layout()
    # Save to file Mesh of edge
    # plt.savefig(f'function_interation_{it}.eps', format='eps', dpi=1200)

    if Show:
        plt.show()