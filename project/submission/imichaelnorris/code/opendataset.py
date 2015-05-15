import downloaddata
import networkx
from rdp import RDPG
from blocks import make_graph_from_blocks
import time
#graph = networkx.read_graphml(open('data/c.elegans_neural.male_1.graphml'))

#parsed from worm atlas http://www.wormatlas.org/neuronalwiring.html
#graph = networkx.read_adjlist('data/wormatlas.graphml')

graph = networkx.read_adjlist('data/wormatlas_downloaded.graphml')
print 'opened dataset'

graph = make_graph_from_blocks(graph)

def clustering_coefficient(graph):
    '''Pavlovic et al. Stochastic Blockmodeling of the Moduels and Core of
    the C. elegans Connectome.  Page 4'''
    numerator = 0.0
    #print 'numerator'
    begin = time.time()
    '''
    for i in sorted(graph.edge):
        for j in sorted(graph.edge):
            for k in sorted(graph.edge):
                if j in graph.edge[i] and \
                   k in graph.edge[j] and \
                   k in graph.edge[i]:
                   numerator += 1
    '''
    for start in graph.edge:
        for mid in graph.edge[start]:
            for end in graph.edge[mid]:
                if start in graph.edge[end]:
                    numerator += 1
    denominator = 0.0

    end = time.time()
    #print 'numerator took', end - begin
    begin = end
    '''
    for i in sorted(graph.edge):
        for k in sorted(graph.edge):
            if k == i:
                continue
            #for j in range(k+1,len(graph.nodes())):
            for j in sorted(graph.edge):
                if j <= k:
                    continue
                if j == i:
                    continue

                if j in graph.edge[i] and k in graph.edge[i]:
                    denominator += 1
    '''
    for node in sorted(graph.edge):
        denominator += len(graph.edge[node])
    end = time.time() 
    #print "denominator took", end-begin
    #print 'num, dem', numerator, denominator
    denominator *= 2.0
    if denominator == 0:
        return 100000
    Cn = numerator / denominator
    return Cn

#print clustering_coefficient(graph)
#if __name__ == '__main__':
    #print clustering_coefficient(graph2)
