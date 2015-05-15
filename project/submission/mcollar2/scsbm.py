#==============================================================================#
# SC-SBM
# ...
#==============================================================================#

# Imports

from numpy import *

import scipy
import scipy.linalg
import scipy.misc
import scipy.io

#==============================================================================#

def sc_ase( Adj, d=2 ):
    
    (U, sVec, tmp) = scipy.linalg.svd( Adj )
    Ud = U[:, :d]
    S = diag( sVec )
    Sd = S[:d, :d]

    Sroot = scipy.linalg.sqrtm( Sd )
    Xhat = dot( Ud, Sroot )
    
    return Xhat

#==============================================================================#

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