'''
Michael Norris
Take C. Elegans connectome and Random Dot Product Graph Mixture Model and try to
find optimal parameters to minimize the difference between graphs generated from
the model and the truth (using various graph metrics, like the clustering
coefficient).
'''

from rdp import RDPGMM
from blocks import blocks
from opendataset import graph as graph_truth, clustering_coefficient
import networkx
import scipy.optimize
import numpy as np
from display_adj_matrix import draw_adj

block_sizes = [len(blocks[block]) for block in sorted(blocks)]
#block_sizes = [4,4,4,4]

def triplet_count(graph):
    triplets = 0

    edges = sorted(graph.edges())

    for i, e1 in enumerate(edges):
        for j, e2 in enumerate(edges[i+1:]):
            for k, e3 in enumerate(edges[i+j+1:]):
                if e1[1] == e2[0] and e2[1] == e3[0] and e3[1] == e1[0]:
                    triplets += 1
    return triplets

def truth_degree_mean_and_stdev(truth):
    '''calculate mean and stdev of node degrees of each block'''
    x = []
    for block in sorted(blocks):
        temp_block = []
        for node in blocks[block]:
            #x.append(len(truth.edge[node]))
            temp_block.append(len(truth.edge[node]))
        x.append( (np.mean(temp_block), np.std(temp_block)))

    return x
    
def block_degree_mean_and_stdev(graph):
    '''gets the stdev of the graph.  If instance of RDPGMM, then we use the adj
    matrix and if it's the truth (networkx) graph, then we will use the blocks
    dict to figure out connectivity of blocks.
    graph: RDPGMM or networkx Graph
    returns: a list with len(block_sizes) tuples with mean and standard 
    deviations of connections.  Each entry is contains the mean and standard
    deviation of the connectivity of each node in that corresponding block.'''
    if type(graph) == RDPGMM:
        x = []
        i= 0
        for block_size in block_sizes:
            temp_block = []
            for node in blocks[block]:
                temp_block.append(sum(graph.adj[i]))

            #x.append( (np.mean(temp_block), np.std(temp_block)))
            x.append(np.mean(temp_block))
            i += block_size
        return x
    else:
        raise Exception('optimization function given wrong type of argument. \
expected RDPGMM, but got {}'.format(type(graph)))

def block_degrees(graph):
    '''calculates the degrees of blocks.  Used to force blocks to have a certain
amount of connections.'''
    if type(graph) == RDPGMM:
        x = []
        i = 0
        for block_size in block_sizes:
            temp_block = []
            x.append(sum(sum(graph.adj[i:i+block_size,:]))/block_size)

            i += block_size
        return x
    else:
        raise Exception('block_degrees() given wrong type of argument. \
expected RDPGMM, but got {}'.format(type(graph)))


def mean_stdev_metric(list1, list2):
    '''Calculates the metric comparing the mean and standard deviation of node
    degrees of partitioned blocks of the RDP-MM graphs.
    returns sum of absolute value of difference between means, where the
    differences of means are scaled by block size,
                and the absolute value of differences in standard deviation,
                where standard deviation is scaled by 
    '''
    m = 0.0


def node_degree(graph):
    '''degrees of all nodes in the graph'''
    return [len(graph.edge[i]) for i in graph.edge]

#alphas = [.6, 1.0, 1.0, 1.0]
#alphas = [2.0] * len(block_sizes)
#.8 to 1.3??
#1.13, 1.40

#variable
#alphas = [1.25]*len(block_sizes)
alphas = [1.25, 1.25, 1.25, 1.5, 0.8,  1.25, 1.25, 1.25]
shape = (sum(block_sizes), sum(block_sizes))

metric = clustering_coefficient
#metric = block_degree_mean_and_stdev
#minimize = lambda x, m=metric: m(RDPGMM(shape, x, block_sizes))- m(graph_truth)

last_value = [0]

#to memoize the graph_truth under a function
truth_memoize_under_function = {}
def minimize(x, m):
    '''function to minimize the metric m on the RDPGMM graph with alphas=x.'''

    for i in x:
        if i < 0:
            return 100000
        if i > 3.0:
            return 100000
    if not m in truth_memoize_under_function:
        #truth_memoize_under_function[m] = m(graph_truth)
        truth_memoize_under_function[m] = m(graph_truth)#\
        truth_memoize_under_function[truth_degree_mean_and_stdev] = \
                                    truth_degree_mean_and_stdev(graph_truth)
            #truth_degree_mean_and_stdev(graph_truth)
    _truth_value = truth_memoize_under_function[m]
    _truth_mean_and_stdev = \
                     truth_memoize_under_function[truth_degree_mean_and_stdev]
    print "args", x
    #mean_diff_weight = .25
    #std_diff_weight = .3
    #value = abs(sum(np.array(m(RDPGMM(shape, x, block_sizes))) - \
    #                np.array(_truth_value)))
    #value = value[0] * mean_diff_weight + value[1] * std_diff_weight
    #print 'evaluation', value
    the_graph = RDPGMM(shape, x, block_sizes)
    last_value[0] = m(the_graph.graph)
    value = abs(last_value[0] - _truth_value)
    
    #number of nodes in each block
    #block_degs = block_degrees(the_graph)
    block_degs = block_degree_mean_and_stdev(the_graph)
    #block_degs = [mean for mean, std in block_degs]

    #make sure that each block's total number of edge degrees are within one
    #standard deviation of the mean of the original graph.
    #this may help to prevent the graph from getting out of hand.

    for block_x, _truth_mean_std in zip(block_degs, _truth_mean_and_stdev):
        _truth_mean, _truth_std = _truth_mean_std
        print block_x, _truth_std, _truth_mean
        #if it's outside of the range, penalize
        if block_x >= _truth_std + _truth_mean:
            value += block_x - (_truth_std + _truth_mean)
        elif block_x <= _truth_std + _truth_mean:
            value += _truth_std + _truth_mean - block_x
        #if (block_x <= _truth_std + _truth_mean) and \
        #   (block_x >= _truth_std - _truth_mean):
        #    pass
        #else:


    print 'val', value
    return value

#optimizing

min_cc = lambda x: minimize(x, metric)

optimize_function = scipy.optimize.fmin_cg

#best
'''
truth metric value {<function truth_degree_mean_and_stdev at 0x3bea6e0>:
[(12.800000000000001, 5.2211109928826449), (25.65625, 7.5771423331952796),
(8.8125, 3.8603553916705646), (48.333333333333336, 9.1772665986241364), (73.5,
15.945218719101975), (22.115384615384617, 5.7333499609868515),
(13.36986301369863, 5.4105600685333259), (18.199999999999999,
2.4549270186029295)], <function clustering_coefficient at 0x3a41e60>:
2.805596005824839}

0.97330373  1.0891615   1.06504759  1.12435104  0.80715165  1.24501222
  1.14101498  0.83222614

'''

#optimize_function = scipy.optimize.fmin

#3.58665959  14.86751322   6.83829878   4.39188916   3.27488678
#0.98334403   1.00087729   0.95337238
#optimize_function = scipy.optimize.fmin_powell

result = optimize_function(min_cc, alphas, epsilon=0.05, maxiter=3)
print 'done, best match', result
print 'truth metric value', truth_memoize_under_function
print 'RDP graph metric value', last_value[0]

best_match = RDPGMM(shape, result, block_sizes)

draw_adj(best_match)

#temp_graph = RDPGMM(shape, alphas, block_sizes)
#print 'temp edges:', len(temp_graph.graph.edges())
#print 'true edges:', len(graph_truth.edges())
#print temp_graph.graph.adjacency_list()
#print 'temp cc', clustering_coefficient(temp_graph.graph)

#print 'temp triplets', triplet_count(temp_graph.graph)
#print 'truth triplets', triplet_count(graph_truth)

#networkx.draw(temp_graph.graph)
#print 'truth cc', clustering_coefficient(graph_truth)
