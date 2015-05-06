
# coding: utf-8

# #StatConn HW1
# Find where $k$-means on the adjacency matrix doesn't "work".
# 
# One reasonable approach at defining "working" is to generate a sample from an SBM, and then see if the labels obtained by $k$-means agree with the "true" block membership via ARI. We'll assume for simplicity that $k$ is known in advance.

# In[48]:

# Imports

from numpy import *

import scipy
import scipy.linalg
import scipy.misc
import scipy.io
import scipy.stats

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

from sklearn import metrics
from sklearn.cluster import KMeans

import time

get_ipython().magic('pylab inline')
pylab.rcParams['savefig.dpi'] = 150


# ##Generate samples from an SBM

# In[51]:

def sc_sbm( rho, P, n ):
    '''sc_sbm : Generate a sample from a (undirected, no-self-loop) SBM.
            rho - Block membership probability vector
            P - Block connection probability matrix
            n - Number of vertices in graph'''
    
    # -- Setup
    k = rho.shape[0]
    membership = zeros( (0) )
    adjacency = zeros( (n, n) )
    
    # -- Generate block membership
    # Make sure rho is normalized
    rho = (1 / sum(rho)) * rho
    # Generate using multinomial
    mems_agg = random.multinomial( n, rho )
    for i in arange( k ):
        for j in arange( mems_agg[i] ):
            membership = r_[ membership, i ]
    
    # -- Generate adjacency
    for i in arange( n ):
        for j in arange( i-1 ):
            u = random.uniform()
            if u < P[ membership[i], membership[j] ]:
                adjacency[ i, j ] = 1
                adjacency[ j, i ] = 1
    
    return ( adjacency, membership )


# In[52]:

# Testing
P = 0.1 * ones( (3,3) ) + 0.3 * eye( 3 )

(A, m) = gen_sbm( array([0.1, 0.3, 0.6]), P, 100 )
plt.set_cmap( 'gray' )
plt.imshow( 1-A, interpolation='nearest' )


# ##Cluster vertices with $k$-means

# In[81]:

# -- Parameters
k = 3
n = 100

rho = random.uniform( size=(k) )
rho = (1 / sum(rho)) * rho

withinP = random.uniform( low=0.2, high=0.6, size=(k) )
P = 0.1 * ones( (k,k) ) + diag( withinP )

nIters = 1000

ariWorst = 1.0
mTrueWorst = []
mEstWorst = []
AWorst = []

aris = []

for i in arange( nIters ):

    # -- Generate
    (A, mTrue) = sc_sbm( rho, P, n )

    # -- Cluster
    km = KMeans( n_clusters=k )
    mEst = km.fit_predict( A )

    # -- Validate
    ariCur = metrics.adjusted_rand_score( mTrue, mEst )
    aris.append( ariCur )
    if ariCur < ariWorst:
        ariWorst = ariCur
        mTrueWorst = mTrue
        mEstWorst = mEst
        AWorst = A
    
    if mod( i, 10 ) == 0:
        print( i, end='... ' )

print('Done!')


# In[100]:

# Display results
print( 'Worst ARI:', ariWorst )

plt.hist( array( aris ), bins=30 )
plt.xlabel( 'ARI' )
plt.ylabel( 'Count' )
plt.title( "Performance distribution for these params." )
plt.show()

plt.figure( figsize=(6,3) )
plt.plot( arange( n ), mTrueWorst, 'k+', markersize=10)
plt.plot( arange( n ), mEstWorst, 'rx', markersize=10)
(x1,x2,y1,y2) = plt.axis()
plt.axis( (x1,x2,-0.5,k-0.5) )
plt.xlabel( 'Vertex' )
plt.ylabel( 'Class' )
plt.legend( ('True', 'Estimated'), loc='upper left' )
plt.title( 'Labeling for worst performance' )
plt.show()

plt.set_cmap( 'gray' )
plt.imshow( 1 - AWorst, interpolation='nearest' )
plt.title( 'Adjacency' )
plt.title( 'Adjacency for worst performance' )
plt.show()


# Note: The figures above were generated with the block connection matrix
# 
# $$ P = \left( \begin{array}{ccc}
# 0.551 & 0.1 & 0.1 \\
# 0.1 & 0.497 & 0.1 \\
# 0.1 & 0.1 & 0.303 \end{array} \right)$$

# ##Debugging miscellanea

# In[49]:

whos


# In[89]:

P

