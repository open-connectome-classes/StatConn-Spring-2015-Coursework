#==============================================================================#
# STATCONN: FINAL PROJECT
# Max Collard
#==============================================================================#

# == IMPORTS

import xml.etree.ElementTree as ET

from numpy import *

import scipy
import scipy.linalg
import scipy.io
import scipy.integrate
import scipy.optimize

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

from sklearn import metrics
from sklearn.cluster import KMeans

import time
import os
import sys

# Custom
import scsbm
import scba
import scdata

#==============================================================================#

print( "STATCONN FINAL -- MAX COLLARD" )
print( "=============================" )
print( "HEAD'S UP: This might take a while. To give you an idea:" )
print( "    SBM     : Progress runs from 20 to 250, increments of 10." )
print( "    BA      : Ditto." )
print( "    Worm+BA : n runs from 0 to 250, increments of 10; each has progress run from 0 to 100 (shows 80)." )
print()
sys.stdout.flush()

fig = plt.figure()

## == SBM

print( 'SBM', end=': ' )
sys.stdout.flush()

# -- Parameters

# Determined
k = 4
ns = arange( 20, 251, 10 )
nn = ns.shape[0]
nIters = 200

# Random
rho = random.uniform( size=(k) )
rho = (1 / sum(rho)) * rho

withinP = random.uniform( low=0.2, high=0.6, size=(k) )
P = 0.1 * ones( (k,k) ) + diag( withinP )

aris = zeros( (nn, nIters) )

# -- Loop

for curNIdx in arange( nn ):

    n = ns[curNIdx]
    #print('n=', n, end=': ')
    print(n, end='...')
    sys.stdout.flush()
    
    for curIter in arange( nIters ):

        # -- Generate
        (A, labelTrue) = scsbm.sc_sbm( rho, P, n )

        # -- Estimate
        Xhat = scsbm.sc_ase( A, k )

        # -- Cluster
        km = KMeans( n_clusters=k )
        labelEst = km.fit_predict( Xhat )

        # -- Validate
        ariCur = metrics.adjusted_rand_score( labelTrue, labelEst )
        aris[curNIdx, curIter] = ariCur

        #if mod( curIter, 10 ) == 0:
        #    print( curIter, end='... ' )
        
    #print()

print('Done!')
print()
sys.stdout.flush()

# -- Plot

# Dummy check for figure directory
figDir = 'Figures'
if not os.path.exists( figDir ):
    os.makedirs( figDir )
    
nBinsHist = 50

histAll = zeros( (nn, nBinsHist) )

for curNIdx in arange( nn ):
    (curHist, tmp1, tmp2) = plt.hist( abs(aris[curNIdx, :]), bins=nBinsHist, range=array([0,1]) )
    histAll[curNIdx, :] = (1 / nIters) * curHist

plt.imshow( flipud( histAll.T ), interpolation='nearest', extent=(ns[0]-5, ns[-1]+5, 0, 1), aspect='auto' )
cb = plt.colorbar()
plt.title('ASE-MSE performance on SBM')
plt.xlabel('# of vertices')
plt.ylabel('ARI')
cb.set_label('Fraction of trials')

fn = 'Figures/sbm-perf-{0}.pdf'.format( random.randint(0, high=1000) )
plt.savefig( fn )
plt.show()

fig.clf()

#==============================================================================#

# == BA

print( 'BA', end=': ' )
sys.stdout.flush()

# -- Parameters

k = 4
Nstart = 5
Nk = Nstart * ones( (k) )

nConns = 3

nsBA = arange( 20, 251, 10 )
nnBA = nsBA.shape[0]
nItersBA = 200

arisBA = zeros( (nnBA, nItersBA) )

# -- Loop

for curNIdx in arange( nnBA ):

    n = nsBA[curNIdx]
    print(n, end='...')
    sys.stdout.flush()
    #print('n=', n, end=': ')
    
    for curIter in arange( nItersBA ):

        # -- Generate
        (A, labelTrue) = scba.start_graph( Nk )
        (A, labelTrue) = scba.ba_iterate( A, labelTrue, iters = n - sum( Nk ), conns = nConns, debug = False )

        # -- Estimate
        Xhat = scsbm.sc_ase( A, k )

        # -- Cluster
        km = KMeans( n_clusters=k )
        labelEst = km.fit_predict( Xhat )

        # -- Validate
        ariCur = metrics.adjusted_rand_score( labelTrue, labelEst )
        arisBA[curNIdx, curIter] = ariCur

        #if mod( curIter, 20 ) == 0:
        #    print( curIter, end='... ' )
        
    #print()

print('Done!')
print()
sys.stdout.flush()

# -- Plotting

nBinsHistBA = 50

histAllBA = zeros( (nnBA, nBinsHistBA) )

for curNIdx in arange( nnBA ):
    (curHist, tmp1, tmp2) = plt.hist( abs(arisBA[curNIdx, :]), bins=nBinsHistBA, range=array([0,1]) )
    histAllBA[curNIdx, :] = (1 / nItersBA) * curHist

plt.imshow( flipud( histAllBA.T ), interpolation='nearest',
            extent=(nsBA[0]-5, nsBA[-1]+5, 0, 1), aspect='auto' )
cb = plt.colorbar()
plt.title('ASE-MSE performance on BA')
plt.xlabel('# of vertices')
plt.ylabel('ARI')
cb.set_label('Fraction of trials')

fn = 'Figures/ba-perf-{0}.pdf'.format( random.randint(0, high=1000) )
plt.savefig( fn )
plt.show()

fig.clf()

#==============================================================================#

# == Worm + BA

print( 'Worm+BA:' )
sys.stdout.flush()

# -- Import GraphML data

#dataUrl = 'https://itsthefedora.github.io/jhu/statconn/c-elegans-male-1.graphml'
dataUrl = 'https://www.dropbox.com/sh/idt3d0gylplyo31/AABy_BCg3yTxfcEDDo0yOCS_a/c-elegans-male-1.graphml?dl=1'
graphFN = 'worm.graphml'
scdata.sc_download( dataUrl, graphFN )

tree = ET.parse( graphFN )
root = tree.getroot()
graphRoot = root.findall( 'graph' )[0]

Aworm = scdata.sc_graphml_adj( graphRoot, symmetrize = True )

# -- Parameters

k = 4

nsWorm = arange( 0, 251, 10 )
nnWorm = nsWorm.shape[0]
nItersWorm = 100

nConns = 3

nVertexWorm = Aworm.shape[0]
arisWorm = zeros( (nnWorm, nItersWorm) )

# -- Determine ground truth

XhatTrue = scsbm.sc_ase( Aworm, k )
km = KMeans( n_clusters=k )
labelTrue = km.fit_predict( XhatTrue )

# -- Loop

for curNIdx in arange( nnWorm ):

    n = nsWorm[curNIdx]
    print('n=', n, end=': ')
    
    for curIter in arange( nItersWorm ):

        # -- Generate
        (A, tmp) = scba.ba_iterate( Aworm, zeros( (nVertexWorm) ),
                                    iters = n, conns = nConns, debug = False )

        # -- Embned
        Xhat = scsbm.sc_ase( A, k )

        # -- Cluster
        km = KMeans( n_clusters=k )
        labelEst = km.fit_predict( Xhat )
        
        # -- Validate
        ariCur = metrics.adjusted_rand_score( labelTrue, labelEst[:nVertexWorm] )
        arisWorm[curNIdx, curIter] = ariCur

        if mod( curIter, 20 ) == 0:
            print( curIter, end='... ' )
            sys.stdout.flush()
        
    print()

print('Done!')
sys.stdout.flush()

# -- Plotting

nBinsHistWorm = 50

histAllWorm = zeros( (nnWorm, nBinsHistWorm) )

for curNIdx in arange( nnWorm ):
    (curHist, tmp1, tmp2) = plt.hist( abs(arisWorm[curNIdx, :]), bins=nBinsHistWorm, range=array([0,1]) )
    histAllWorm[curNIdx, :] = (1 / nItersBA) * curHist

plt.imshow( flipud( histAllWorm.T ), interpolation='nearest',
            extent=(nsWorm[0]-5, nsWorm[-1]+5, 0, 1), aspect='auto' )
cb = plt.colorbar()
plt.title('ASE-MSE performance on Worm+BA')
plt.xlabel('# of added vertices')
plt.ylabel('ARI')
cb.set_label('Fraction of trials')

fn = 'Figures/wormba-perf-{0}.pdf'.format( random.randint(0, high=1000) )
plt.savefig( fn )
plt.show()

fig.clf()


print()
print( 'Have a nice day! :)' )
sys.stdout.flush()






































