'''
Random Dot Product Graph Model
'''
import numpy as np
from scipy.stats import bernoulli
import networkx as nx
import pylab

class RDPG(object):
    def __init__(self, shape, alpha=1.0):
        assert(type(shape) == tuple)
        assert(len(shape) == 2)
        assert(shape[0] == shape[1])
        self.alpha = alpha
        self._gen_adjacency_matrix(shape)
        self._graph = None
        self.shape = shape

    def _gen_adjacency_matrix(self, shape):
        self.adj = np.zeros(shape)
        #directed RDPG has X, Y matrices.  See RDPGM for SN Young & Scheiemerman
        #pg 3
        #uniform by pg 2
        low = 0.0
        high = 1.0
        d = shape[0]
        self.shape = shape
        
        #set random seed for reproducability
        self.X = np.random.uniform(low=low, high=high, size=shape)
        self.Y = np.random.uniform(low=low, high=high, size=shape)

        self.X *= 1.0 / np.sqrt(d)
        self.Y *= 1.0 / np.sqrt(d)

        self.X **= self.alpha
        self.Y **= self.alpha

        while not self._assert_dot_products_are_probabilities(self.X, self.Y):
            pass

        for i in range(shape[0]):
            for j in range(shape[1]):
                if i == j:
                    continue

                self.adj[i][j] = bernoulli.rvs(self.X[i].dot(self.Y[j]))

    def _assert_dot_products_are_probabilities(self, mat1, mat2):
        for i,col1 in enumerate(mat1.transpose()):
            for j,col2 in enumerate(mat2.transpose()):
                #if we violate this, we need to regenerate
                if i == j:
                    continue

                if ( abs(col1.dot(col2)) > 1.0):
                    print "dot product <= 1.0 condition violated: regenerating"
                    return False
        return True

    @property
    def graph(self):
        if self._graph == None:
            self._graph = nx.Graph(self.adj)
        return self._graph

    def embed(self):
        #also see sklearn.manifold.SpectralEmbedding
        U, S, V = np.linalg.svd(self.adj)
        #what happens when n != d
        X = U.dot(np.sqrt(S.transpose() * np.eye(self.shape[0])))
        Y = V.dot(np.sqrt(S.transpose() * np.eye(self.shape[0])))
        return X, Y

    def clustering_coefficient(self):
        '''Pavlovic et al. Stochastic Blockmodeling of the Moduels and Core of
        the C. elegans Connectome.  Page 4'''
        graph = self.graph
        numerator = 0.0
        for i in sorted(graph.edge):
            for j in sorted(graph.edge):
                for k in sorted(graph.edge):
                    #numerator += graph.adj[i][j] * \
                    #             graph.adj[j][k] * \
                    #             graph.adj[i][k]
                    if j in graph.edge[i] and \
                        k in graph.edge[j] and \
                        k in graph.edge[i]:
                        numerator += 1
        denominator = 0.0

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
                #denominator += graph.adj[i][j] * graph.adj[i][k]
                if j in graph.edge[i] and k in graph.edge[i]:
                    denominator += 1
    
        denominator *= 2.0
        Cn = numerator / denominator
        return Cn

    #This is a CC directly from the adj matrix but we can get the networkx graph
    #representation first, which I prefer.
    def __clustering_coefficient(self):
        '''Pavlovic et al. Stochastic Blockmodeling of the Moduels and Core of
        the C. elegans Connectome.  Page 4'''
        numerator = 0.0
        for i in range(self.shape[0]):
            for j in range(self.shape[0]):
                for k in range(self.shape[0]):
                    numerator += self.adj[i][j] * \
                                 self.adj[j][k] * \
                                 self.adj[i][k]
        denominator = 0.0

        for i in range(self.shape[0]):
            for k in range(self.shape[0]):
                if k == i:
                    continue
                for j in range(k+1,self.shape[0]):
                    if j == i:
                        continue
                denominator += self.adj[i][j] * self.adj[i][k]
        
        denominator *= 2.0
        Cn = numerator / denominator
        return Cn

class RDPGMM(RDPG):
    def __init__(self, shape, alpha=1.0, block_sizes=None):
        assert(type(shape) == tuple)
        assert(len(shape) == 2)
        assert(shape[0] == shape[1])
        self.alpha = alpha
        self.block_sizes = block_sizes
        self._gen_adjacency_matrix(shape, block_sizes)
        self._graph = None
        self.shape = shape

    def _gen_adjacency_matrix(self, shape, block_sizes=None):
        '''Generate RDP of given shape with list of block sizes
        shape: 2-tuple, shape[0] == shape[1] 
        block_sizes: list of block sizes [b1, b2, b3], 
                     sum(block_sizes) == shape[0]'''

        #if no list
        if not block_sizes:
            block_sizes = [shape[0]]

        self.adj = np.zeros(shape)
        #directed RDPG has X, Y matrices.  See RDPGM for SN Young & Scheiemerman
        #pg 3
        #uniform by pg 2
        low = 0.0
        high = 1.0
        d = shape[0]
        self.shape = shape

        self.X = np.zeros(shape)
        self.Y = np.zeros(shape)

        #self.X = np.random.uniform(low=low, high=high, size=shape)
        #self.Y = np.random.uniform(low=low, high=high, size=shape)
        x = 0
        y = 0

        assert sum(block_sizes) == self.shape[0]
        for i,block_size in enumerate(block_sizes):
            _shape = (block_size, d)
            self.X[x:x+block_size] = np.random.uniform(low=low,
                                                         high=high,size=_shape)
            #self.Y[y:y+block_size] = np.random.uniform(low=low,
            #                                             high=high,size=_shape)

            self.X[x:x+block_size] *= 1.0 / np.sqrt(d)
            #self.Y[y:y+block_size] *= 1.0 / np.sqrt(d)

            self.X[x:x+block_size] **= self.alpha[i]
            #self.Y[y:y+block_size] **= self.alpha[i]

            x += block_size
            y += block_size

        self.Y = self.X

        assert x == self.shape[0]
        assert y == self.shape[1]

        retries = 10

        #i = 0
        #while not self._assert_dot_products_are_probabilities(self.X, self.Y)\
        #     and i < retries:
        #     i += 1

        #if i >= retries:
        #    print "warning: paramater caused dot product to be greater than 1"

        #for i in range(shape[0]):
        #    for j in range(shape[1]):
        #        if i == j:
        #            continue
        #        self.adj[i][j] = bernoulli.rvs(min(self.X.transpose()[i].dot(
        #                                       self.Y.transpose()[j]),1))
        self.adj = self.X.dot(self.Y)
        self.adj = np.clip(self.adj, a_min=0, a_max=1.0)
        self.adj = bernoulli.rvs(self.adj)

if __name__ == "__main__":
    #graph = RDPG((16, 16))
    graph = RDPGMM((16, 16), alpha=[1.0, 1.1, 1.2, 1.3], block_sizes=[4,4,4,4])
    #print graph.clustering_coefficient()
    #print(type(graph.graph))

    #print graph.embed()
    #nx.draw(graph.graph)
    #pylab.show()

