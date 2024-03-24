# Copyright (c) 2018 Gdansk University of Technology
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# Author: Piotr Kowalczyk
# Project homepage: https://github.com/PioKow/GRPF
#
# Python implementation: Michał Wiliński


import numpy as np
from scipy.spatial import Delaunay
import triangle as tr
import matplotlib.pyplot as plt

import matplotlib.pyplot

from vinq import vinq
from rect_dom import rect_dom
from find_next_node import find_next_node
from disk_dom import disk_dom
from vis import vis
from candidate_edges_Q import find_candidate_edges


def grpf():

    from analysis_parameters import ItMax, NodesMax, Tol, SkinnyTriangle, visual, NewNodesCoord
    from fun import fun
        
    # Initialize variables
    NodesCoord = np.array([]).reshape(0, 2)
    NrOfNodes = 0
    FuntionValues = np.zeros((ItMax * NodesMax), dtype=complex)
    Quadrants = np.zeros((0, 1))

    # Loop that analyzes the mesh, finds 
    it = 0
    while it < ItMax and NrOfNodes < NodesMax:
        it += 1
    
        NodesCoord = np.vstack([NodesCoord, NewNodesCoord])

        # Get node indicies for quantum GRPF (order by x and y)
        dtype = [('x', float), ('y', float)]
        tempNodesCoord = np.core.records.fromarrays([NodesCoord[:,0], NodesCoord[:,1]], names='x, y', formats='f, f')
        node_indicies = np.argsort(tempNodesCoord, order=['x', 'y']).flatten()
        
        print(f'Evaluation of the function in {len(NewNodesCoord)} new points...')
    
        # Evaluation of quadrants of function values in new coordinates 
        Quadrants = np.vstack([Quadrants, np.zeros((len(NewNodesCoord), 1))])
        for Node in range(NrOfNodes, NrOfNodes + len(NewNodesCoord)):
            z = NodesCoord[Node, 0] + 1j * NodesCoord[Node, 1]
            FuntionValues[Node] = fun(z)
            Quadrants[Node] = vinq(FuntionValues[Node])
    
        NewNodesCoord = np.empty((0, 2))
        NrOfNodes = len(NodesCoord)

        print(f'Triangulation and analysis of {NrOfNodes} nodes...')

        #Triangulation 
        opts = None
        DT = Delaunay(NodesCoord)
        Elements = DT.simplices  # gets indices of NodesCoord([[x,y], [x,y], ...]) that make triangles([[a,b,c], [a,b,c]]) ->
        # -> NodesCoord[Elements] = [[[x_a, y_a], [x_b, y_b], [x_c, y_c]], [x_a, y_a], [x_b, y_b], [x_c, y_c]], ...]

        NrOfElements = Elements.shape[0]

        # Creation of a list of unique edges from triangles
        Edges = []

        for triangle in Elements:
            Edges.extend([np.sort((triangle[0], triangle[1])), np.sort((triangle[1], triangle[2])), np.sort((triangle[2], triangle[0]))])
        Edges = np.unique(Edges, axis=0)  # indicies of points that make edges ([a, b], [a, b], ...]) -> 

        #Finding candidate edges
        PhasesDiff = np.mod(Quadrants[Edges[:, 0]] - Quadrants[Edges[:, 1]], 4).flatten()

        CandidateEdges = []
        
        if it == 1:
            # Quantum part

            max_tries = 10
            are_edges_true = False
            for tries in range(max_tries):
                print(f'tries {tries}')
                
                # find grid width (counted in number of points in the first line in quantum ordering -> leftmost straight line)
                grid_width = 0
                for i in range(len(node_indicies)):
                    if i != 0 and NodesCoord[node_indicies[0], 0] != NodesCoord[node_indicies[i], 0]:
                        grid_width = i
                        break
        
                quadrant_arr = Quadrants[:,0].astype(int)
                candidate_edges_Qorder = find_candidate_edges(quadrant_arr[node_indicies], grid_width)  # candidate_edges in order of the quantum algorithm of GRPF
                
                candidate_edges_Corder = np.empty([len(candidate_edges_Qorder), 2], dtype=int)  # candidate_edges in order of the classical algorithm of GRPF
                for i in range(len(candidate_edges_Qorder)):
                    candidate_edges_Corder[i, 0] = node_indicies[candidate_edges_Qorder[i, 0]]
                    candidate_edges_Corder[i, 1] = node_indicies[candidate_edges_Qorder[i, 1]]
    
                are_edges_true = True
                # check if the candidates edges are true
                for Edge in candidate_edges_Corder:
                    if np.absolute(Quadrants[Edge[0]] - Quadrants[Edge[1]]) != 2:
                        are_edges_true = False
                        break
    
                print(f'Edges found by quantum algorithm: {candidate_edges_Corder}')
                print(f'Edges are {are_edges_true}!')
    
                if are_edges_true == True:
                    CandidateEdges = candidate_edges_Corder
                    break

        else:
            CandidateEdges = Edges[(PhasesDiff == 2) | np.isnan(PhasesDiff)]
        
        if len(CandidateEdges) == 0:
            print('No roots in the domain!')
            break
    
        Nodes1OfCandidateEdges = CandidateEdges[:, 0]
        Nodes2OfCandidateEdges = CandidateEdges[:, 1]
    
        CoordinatesOfNodes1OfCandidateEdges = NodesCoord[Nodes1OfCandidateEdges, :]
        CoordinatesOfNodes2OfCandidateEdges = NodesCoord[Nodes2OfCandidateEdges, :]

        # Check if desired accuracy has already been achieved
        CandidateEdgesLengths = np.sqrt(np.sum((CoordinatesOfNodes2OfCandidateEdges - CoordinatesOfNodes1OfCandidateEdges) ** 2, axis=1))
        MinCandidateEdgesLengths = np.min(CandidateEdgesLengths)
        MaxCandidateEdgesLengths = np.max(CandidateEdgesLengths)
        print(f'Candidate edges length min: {MinCandidateEdgesLengths} max: {MaxCandidateEdgesLengths}')
    
        if MaxCandidateEdgesLengths < Tol:
            print(f'Assumed accuracy is achieved in iteration: {it}')
            break

        # We only analyse the edges longer than our length threshold
        Temp = CandidateEdgesLengths > Tol
        ReduCandidateEdges = CandidateEdges[Temp, :]
    
        Temp = np.zeros(NrOfNodes)
        Temp[ReduCandidateEdges[:, 0]] = 1
        Temp[ReduCandidateEdges[:, 1]] = 1
        CandidateNodes = np.asarray(Temp == 1).nonzero()
        CandidateNodes = CandidateNodes[0]
        

        # Find elements that neighbour with Candidate nodes
        ArrayOfCandidateElements = [[i for i, simplex in enumerate(Elements) if vertex in simplex] for vertex in CandidateNodes]
        
        Temp = np.zeros(NrOfElements)
        for part in ArrayOfCandidateElements:
            Temp[part] += 1
        
        # Divide them into groups depending on whether they have a candidate edge, or just a candidate node
        IDOfFirsrZoneCandidateElements = np.asarray(Temp > 1).nonzero()[0]
        IDOfSecondZoneCandidateElements = np.asarray(Temp == 1).nonzero()[0]

        
        NoOfFirsrZoneCandidateElements = len(IDOfFirsrZoneCandidateElements)
        FirsrZoneCandidateElements = Elements[IDOfFirsrZoneCandidateElements, :]
        TempExtraEdges = np.empty((0, 2), dtype=int)

        
        for k in range(NoOfFirsrZoneCandidateElements):
            TempExtraEdges = np.vstack((TempExtraEdges, [FirsrZoneCandidateElements[k, 0], FirsrZoneCandidateElements[k, 1]]))
            TempExtraEdges = np.vstack((TempExtraEdges, [FirsrZoneCandidateElements[k, 1], FirsrZoneCandidateElements[k, 2]]))
            TempExtraEdges = np.vstack((TempExtraEdges, [FirsrZoneCandidateElements[k, 2], FirsrZoneCandidateElements[k, 0]]))

        # Divide the candidate elements by creating a new node in the middle of each edge
        # adding them to the new coordinates list (without repetitions)
        NewNodesCoord = np.sum(NodesCoord[TempExtraEdges[0, :], :], axis=0) / 2
        NewNodesCoord = NewNodesCoord.reshape(1, 2)
        for k in range(1, 3 * NoOfFirsrZoneCandidateElements):
            CoordOfTempEdgeNode1 = NodesCoord[TempExtraEdges[k, 0], :]
            CoordOfTempEdgeNode2 = NodesCoord[TempExtraEdges[k, 1], :]
            TempNodeCoord = (CoordOfTempEdgeNode1 + CoordOfTempEdgeNode2) / 2
            TempEdgeLength = np.sqrt(np.sum((CoordOfTempEdgeNode2 - CoordOfTempEdgeNode1) ** 2))

            if TempEdgeLength > Tol:
                DistNodes = np.sqrt((NewNodesCoord[:, 0] - TempNodeCoord[0]) ** 2 + (NewNodesCoord[:, 1] - TempNodeCoord[1]) ** 2)
                if np.sum(DistNodes < 2 * np.finfo(float).eps) == 0:
                    NewNodesCoord = np.vstack((NewNodesCoord, TempNodeCoord))
                    

        
        # removing the first new node if the edge is too short
        CoordOfTempEdgeNode1 = NodesCoord[TempExtraEdges[0, 0], :]
        CoordOfTempEdgeNode2 = NodesCoord[TempExtraEdges[0, 1], :]
        TempEdgeLength = np.sqrt(np.sum((CoordOfTempEdgeNode2 - CoordOfTempEdgeNode1) ** 2))
        if TempEdgeLength < Tol:
            NewNodesCoord = np.delete(NewNodesCoord, 0, axis=0)
    
        NoOfSecondZoneCandidateElements = len(IDOfSecondZoneCandidateElements)
        SecondZoneCandidateElements = Elements[IDOfSecondZoneCandidateElements, :]


        # Adding singular nodes for the second zone candidate elements,
        # and only if it wouldn't make the triangles "too skinny"
        for k in range(NoOfSecondZoneCandidateElements):
            NodesInTempElement = SecondZoneCandidateElements[k, :]
            Node1Coord = NodesCoord[NodesInTempElement[0], :]
            Node2Coord = NodesCoord[NodesInTempElement[1], :]
            Node3Coord = NodesCoord[NodesInTempElement[2], :]
            TempLengths = np.zeros(3)
            TempLengths[0] = np.sqrt(np.sum((Node2Coord - Node1Coord) ** 2))
            TempLengths[1] = np.sqrt(np.sum((Node3Coord - Node2Coord) ** 2))
            TempLengths[2] = np.sqrt(np.sum((Node1Coord - Node3Coord) ** 2))
            if np.max(TempLengths) / np.min(TempLengths) > SkinnyTriangle:
                TempNodeCoord = (Node1Coord + Node2Coord + Node3Coord) / 3
                NewNodesCoord = np.vstack((NewNodesCoord, TempNodeCoord))

        
        if visual == 2:
            vis(NodesCoord, Edges, Quadrants, PhasesDiff)
            pass
        
        print(f'Iteration: {it} done')
        print('----------------------------------------------------------------')

    

    # End of loop
    # Evaluation of regions and verification

    
    # Visualization if visual > 0
    if visual > 0:
        vis(NodesCoord, Edges, Quadrants, PhasesDiff, False)
    
    
    print('Evaluation of regions and verification...')
    
    # Extraction of elements (triangles) that have at least one edge in the set of candidate edges
    ArrayOfCandidateElements = [[i for i, simplex in enumerate(Elements) if (edge[0] in simplex and edge[1] in simplex)] for edge in CandidateEdges]

    Temp = np.zeros(NrOfElements)
    for k in range(len(ArrayOfCandidateElements)):
        Temp[ArrayOfCandidateElements[k]] = 1

    
    IDOfCandidateElements = np.where(Temp == 1)[0]
    NoOfCandidateElements = len(IDOfCandidateElements)
    CandidateElements = Elements[IDOfCandidateElements, :]


    # Extraction of Edges in those candidate elements
    TempEdges = np.zeros((NoOfCandidateElements * 3, 2), dtype=complex)
    
    for k in range(NoOfCandidateElements):
        TempEdges[(k - 1) * 3 + 1, :] = [CandidateElements[k, 0], CandidateElements[k, 1]]
        TempEdges[(k - 1) * 3 + 2, :] = [CandidateElements[k, 1], CandidateElements[k, 2]]
        TempEdges[(k - 1) * 3 + 3, :] = [CandidateElements[k, 2], CandidateElements[k, 0]]

    
    # Reduction of edges to contour
    MultiplicationOfTempEdges = np.zeros(3 * NoOfCandidateElements)
    RevTempEdges = np.fliplr(TempEdges)
    for k in range(3 * NoOfCandidateElements):
        if MultiplicationOfTempEdges[k] == 0:
            NoOfEdge = np.where((RevTempEdges[:, 0] == TempEdges[k, 0]) & (RevTempEdges[:, 1] == TempEdges[k, 1]))[0]
            if len(NoOfEdge) == 0:
                MultiplicationOfTempEdges[k] = 1
            else:
                MultiplicationOfTempEdges[k] = 2
                MultiplicationOfTempEdges[NoOfEdge] = 2
    
    ContourEdges = TempEdges[MultiplicationOfTempEdges == 1, :].astype(int)
    NoOfContourEdges = ContourEdges.shape[0]
    
    # Evaluation of the regions
    NoOfRegions = 1
    Regions = {NoOfRegions: [ContourEdges[0, 0]]}
    RefNode = ContourEdges[0, 1]
    ContourEdges = np.delete(ContourEdges, 0, axis=0)

    #  
    while ContourEdges.shape[0] > 0:
        IndexOfNextEdge = np.where(ContourEdges[:, 0] == RefNode)[0]
        
        if len(IndexOfNextEdge) == 0:
            Regions[NoOfRegions].append(RefNode)
            if ContourEdges.shape[0] > 0:
                NoOfRegions += 1
                Regions[NoOfRegions] = [ContourEdges[0, 0]]
                RefNode = ContourEdges[0, 1]
                ContourEdges = np.delete(ContourEdges, 0, axis=0)
        else:
            if len(IndexOfNextEdge) > 1:
                IndexOfNextEdge
                PrevNode = Regions[NoOfRegions][-1]
                TempNodes = ContourEdges[IndexOfNextEdge, 1].astype(int)
                Index = find_next_node(NodesCoord, PrevNode.astype(int), RefNode.astype(int), TempNodes.astype(int))
                IndexOfNextEdge = IndexOfNextEdge[Index]
            else:
                IndexOfNextEdge = IndexOfNextEdge[0]
                
            NextEdge = ContourEdges[IndexOfNextEdge, :]
            Regions[NoOfRegions].append(ContourEdges[IndexOfNextEdge, 0])
            RefNode = ContourEdges[IndexOfNextEdge, 1]
            ContourEdges = np.delete(ContourEdges, IndexOfNextEdge, axis=0)
    
    # Update the last region
    Regions[NoOfRegions].append(RefNode)


    # Evaluate how many times quadrants change while travelling the QuadrantSequence in the clockwise direction
    # Mean of all the points in the sequence is our approximation of the root/pole
    q = np.empty((NoOfRegions, 1), dtype=int)
    z = np.empty((NoOfRegions, 1), dtype=complex)
    for k in range(1, NoOfRegions+1):
        QuadrantSequence = Quadrants[Regions[k]]

        dQ = QuadrantSequence[1:] - QuadrantSequence[:-1]
        dQ[dQ == 3] = -1
        dQ[dQ == -3] = 1
        dQ[np.abs(dQ) == 2] = np.nan
        q[k-1] = np.sum(dQ) / 4
        z[k-1] = np.mean(NodesCoord[Regions[k], 0]) + 1j * np.mean(NodesCoord[Regions[k], 1])

        print(f'Region: {k} z = {z[k-1]} with q = {q[k-1]}')

        if visual > 0:
            plt.plot(NodesCoord[Regions[k], 0], NodesCoord[Regions[k], 1], linewidth=2, color='c')
    
    if visual > 0:
        plt.show()


    
    # Printing results
    z_root = z[q > 0]
    z_roots_multiplicity = q[q > 0]
    
    NoOfRoots = len(z_root)
    
    print('------------------------------------------------------------------')
    print('Root and its multiplicity: ')
    for i in range(len(z_root)):
        print(z_root[i], "x",  z_roots_multiplicity[i])
    
    z_poles = z[q < 0]
    z_poles_multiplicity = -q[q < 0]
    
    NoOfPoles = len(z_poles)
    
    print('------------------------------------------------------------------')
    print('Poles and its multiplicity: ')
    
    for i in range(len(z_poles)):
        print(z_poles[i], "x" , z_poles_multiplicity[i])

grpf()
