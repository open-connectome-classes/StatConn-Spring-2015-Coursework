#==============================================================================#
# SC-BA
# ...
#==============================================================================#

# Imports

from numpy import *

import scipy
import scipy.linalg
import scipy.io
import scipy.integrate
import scipy.optimize

#==============================================================================#

def start_graph( Nk ):
    '''start_graph DOCSTRING'''
    
    k = Nk.size
    Ntotal = sum( Nk )
    retAdj = zeros( (Ntotal, Ntotal) )
    retLabel = zeros( (Ntotal) )
    
    clusterStarts = r_[0, cumsum(Nk)]
    for iCluster in arange( k ):
        clusterStart = clusterStarts[iCluster]
        for iFrom in arange( clusterStart, clusterStart + Nk[iCluster] ):
            for iTo in arange( clusterStart, clusterStart + Nk[iCluster] ):
                if not (iFrom == iTo):
                    retAdj[iFrom, iTo] = 1.0
            retLabel[iFrom] = iCluster
    
    return (retAdj, retLabel)

#==============================================================================#

# Utility
def label_count( label ):
    '''label_count docstring'''
    
    nLabels = max( label ) + 1
    labels = arange( nLabels )   # Integer ID of labels
    counts = zeros( (nLabels) )  # Number of samples with each label
    
    for i in arange( label.size ):
        counts[label[i]] += 1
    
    return (counts, labels)

def ba_iterate( graph, label, iters = 1, conns = 3, debug = False ):
    '''ba_iterate DOCSTRING'''
    
    # Copy existing graph
    Ntotal = label.size + iters
    
    newLabel = zeros( (Ntotal) )
    newGraph = zeros( (Ntotal, Ntotal) )

    newLabel[:label.size] = label
    newGraph[:label.size, :label.size] = graph
    
    # For each iteration ...
    for curIter in arange( iters ):
        
        NExist = label.size + curIter
        
        # Determine connection probability
        degree = zeros( (NExist) )
        for iVertex in arange( NExist ):
            degree[iVertex] = sum( newGraph[:, iVertex] )
            
        totalDegree = sum( degree )
        pConnect = (1 / totalDegree) * degree
        pCum = r_[0, cumsum(pConnect)]
        
        # Connect
        #randVals = random.rand( conns )
        #willConnect = randVals < pConnect
        
        willConnect = zeros( (NExist) )
        while True:
            pCur = random.rand()
            connectIdx = nonzero( diff( pCum < pCur ) )
            willConnect[connectIdx] = 1
            
            if sum( willConnect ) >= conns:
                break
        
        #print(willConnect)
        
        newGraph[NExist, :NExist] = willConnect
        newGraph[:NExist, NExist] = willConnect
        
        # Choose label
        curLabel = newLabel[:NExist]
        connectLabel = curLabel[nonzero(willConnect)]
        (counts, tmp) = label_count( connectLabel )
        if debug:
            print('Counts: ', counts)

        bestCount = max( counts )
        bestLabel = transpose( nonzero( counts == bestCount ) )
        if debug:
            print('Best labels: ', bestLabel.T)

        newLabel[NExist] = bestLabel[ random.randint( bestLabel.size ) ]
    
    return (newGraph, newLabel)

