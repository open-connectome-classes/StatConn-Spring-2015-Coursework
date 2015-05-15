#==============================================================================#
# SC-DATA
# ...
#==============================================================================#

# Imports

import xml.etree.ElementTree as ET

from numpy import *

import urllib.request, urllib.parse, urllib.error

#==============================================================================#

def sc_download( urlRemote, fnLocal ):
    '''sc_download docstring '''
    
    print('Downloading', urlRemote, end='...')
    urllib.request.urlretrieve(urlRemote, fnLocal)
    print('Done.')

#==============================================================================#
   
def sc_graphml_adj( graphElement, symmetrize = False ):
    '''sc_graphml_adj docstring '''
    
    # Vertices
    nodeDict = dict()
    curNode = 0
    for node in graphElement.findall('node'):
        nodeDict[ node.attrib['id'] ] = curNode
        curNode += 1
    Ntotal = curNode
    
    retAdj = zeros( (Ntotal, Ntotal) )
    
    # Edges
    for edge in graphElement.findall('edge'):
        
        iSource = nodeDict[ edge.attrib['source'] ]
        iTarget = nodeDict[ edge.attrib['target'] ]
        
        retAdj[ iSource, iTarget ] = 1
        if symmetrize:
            retAdj[ iTarget, iSource ] = 1
        
    return retAdj