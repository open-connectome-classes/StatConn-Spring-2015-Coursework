import numpy as np
import networkx as nx
from matplotlib import pyplot, patches
from collections import defaultdict
from scipy import io
import pylab
from blocks import blocks
from rdp import RDPGMM
import matplotlib.pyplot as plt
import matplotlib.cm #colormap

def draw_adjacency_matrix(G, node_order=None, partitions=[], colors=[]):
    """
    - G is a netorkx graph
    - node_order (optional) is a list of nodes, where each node in G
          appears exactly once
    - partitions is a list of node lists, where each node in G appears
          in exactly one node list
    - colors is a list of strings indicating what color each
          partition should be
    If partitions is specified, the same number of colors needs to be
    specified.
    """
    adjacency_matrix = nx.to_numpy_matrix(G, dtype=np.bool, nodelist=node_order)

    #Plot adjacency matrix in toned-down black and white
    fig = pyplot.figure(figsize=(5, 5)) # in inches
    pyplot.imshow(adjacency_matrix,
                  cmap="Greys",
                  interpolation="none")
    
    # The rest is just if you have sorted nodes by a partition and want to
    # highlight the module boundaries
    assert len(partitions) == len(colors)
    ax = pyplot.gca()
    for partition, color in zip(partitions, colors):
        current_idx = 0
        for module in partition:
            ax.add_patch(patches.Rectangle((current_idx, current_idx),
                                          len(module), # Width
                                          len(module), # Height
                                          facecolor="none",
                                          edgecolor=color,
                                          linewidth="1"))
            current_idx += len(module)

def assignmentArray_to_lists(assignment_array):
    by_attribute_value = defaultdict(list)
    for node_index, attribute_value in enumerate(assignment_array):
        by_attribute_value[attribute_value].append(node_index)
    return by_attribute_value.values()

'''
A = io.mmread("adjacency_matrix/Caltech.mtx")
G = nx.from_scipy_sparse_matrix(A)

# Load in array which maps node index to dorm number
# Convert this to a list of lists indicating dorm membership
dorm_assignment = np.genfromtxt("adjacency_matrix/caltech_dorms_blanksInferred.txt", dtype="u4")
dorm_lists = assignmentArray_to_lists(dorm_assignment)

# Create a list of all nodes sorted by dorm, and plot
# adjacency matrix with this ordering
nodes_dorm_ordered = [node for dorm in dorm_lists for node in dorm]
draw_adjacency_matrix(G, nodes_dorm_ordered, [dorm_lists],["blue"])

pylab.show()
'''
def draw_adj(rdp_graph):
    ax = plt.subplot(111, aspect='equal')
    pylab.xlim([0,rdp_graph.shape[0]])
    pylab.ylim([0,rdp_graph.shape[0]])
    ax.imshow(np.ones(rdp_graph.shape) - rdp_graph.adj,
              cmap=matplotlib.cm.Greys_r)
    draw_blocks(rdp_graph.block_sizes, rdp_graph.shape[0], ax)
    plt.show()

def draw_blocks(lines, size, ax):
    longways = np.linspace(0, size, 2)

    i = 0
    for line in lines:
        temp = line
        line = line + i
        i += temp
        
        temp = np.linspace(line,line,2)
        ax.plot(longways, temp, 'blue', temp, longways, 'blue')

if __name__ == '__main__':
    block_sizes = [len(blocks[block]) for block in sorted(blocks)]
    graph = RDPGMM((246,246), [3.1,5.9,10.1,3.5,4.0,7.3,1.2,1.2,.9], block_sizes)
    draw_adj(graph)
    #X,Y = np.meshgrid(np.arange(246), np.arange(246))
    #print graph.adj
    #print graph.adj

    #ax = plt.subplot(111, aspect='equal')
    #ax = plt.subplot(111)
    #ax.imshow(np.ones((246,246)) - graph.adj,cmap = matplotlib.cm.Greys_r)
    #ax.pcolor(X,Y,graph.adj, edgecolor='white', linewidth=2)
    #ax.pcolormesh(X,Y,graph.adj, edgecolor='white', linewidth=2)
    plt.show()
