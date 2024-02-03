# # Converting mechanism from classical to quantum node numbering (OUTDATED)
# find W (width counted in number of points in a snake)
# W = 0 
# for node in range(len(NodesCoord)):
#     if node != 0 and NodesCoord[node, 0] == NodesCoord[0, 0]:
#         W = node
#         break
#
# row = 0  # index of a row
# row_node = 0
# last_row = int(len(NodesCoord) / W) - 1
# #  NodesCoord = k * W + W/2 (k -> number of snake rows + one long row; W/2 -> number of nodes in the half row (bottom nodes))
# for node in range(len(NodesCoord)):
#     #plt.triplot(NodesCoord[:,0], NodesCoord[:,1], Elements)
#     plt.plot(NodesCoord[:,0], NodesCoord[:,1], 'o')
#
#     if row <= last_row:
#         # snake right
#         edgeA = np.stack((NodesCoord[node, :], NodesCoord[node, :]))
#         if row_node + 1 != W and node + 1 < len(NodesCoord):
#             edgeA = np.stack((NodesCoord[node, :], NodesCoord[node + 1, :]))
#         plt.plot(edgeA[:,0] , edgeA[:,1])
#
#     if row < last_row:
#         if row_node % 2 != 0:
#             # single left-up
#             if row_node != 0: 
#            
#                 edgeB = np.stack((NodesCoord[node, :], NodesCoord[node, :]))
#                 if node + W - 1 < len(NodesCoord) and row_node != 0:
#                     edgeB = np.stack((NodesCoord[node, :], NodesCoord[node + W - 1, :]))
#                 plt.plot(edgeB[:,0] , edgeB[:,1])
#
#             # single right-up
#             if row_node != W - 1: 
#                 edgeB = np.stack((NodesCoord[node, :], NodesCoord[node, :]))
#                 if node + W + 1 < len(NodesCoord) and row_node != 0:
#                     edgeB = np.stack((NodesCoord[node, :], NodesCoord[node + W + 1, :]))
#                 plt.plot(edgeB[:,0] , edgeB[:,1])
#
#         # single up
#         edgeB = np.stack((NodesCoord[node, :], NodesCoord[node, :]))
#         if node + W < len(NodesCoord):
#             edgeB = np.stack((NodesCoord[node, :], NodesCoord[node + W, :]))
#         plt.plot(edgeB[:,0] , edgeB[:,1])
#
#     if row == last_row + 1:
#         pass
#     row_node += 1
#     if row_node == W:
#         row += 1
#         row_node = 0
# plt.show()