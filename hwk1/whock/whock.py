from sklearn.cluster import KMeans
import numpy as np
from matplotlib import pyplot as plt
import networkx as nx

def playwithkmeans(n=10,p=0.5):
    '''randomly generate a Erdos-Renyi graph with n nodes and p prob of connection
    clust number is, for now, hard-coded to 3 clusts. Raw graph and graph labeled
    according to kmeans clustering are output to see if kmeans "fails" '''
    G = nx.erdos_renyi_graph(n,p)
    pos = nx.circular_layout(G)
    adjmat = nx.to_numpy_matrix(G)
    km = KMeans(n_clusters=3)
    kmfit = km.fit(adjmat)
    l = kmfit.labels_ 
    c1 = []
    c2 = []
    c3 = []
    for i,x in enumerate(l):
        if x == 0:
            c1.append(i)
        if x == 1:
            c2.append(i)
        if x == 2:
            c3.append(i)
    nx.draw_networkx_nodes(G,pos,
                       nodelist=c1,
                       node_color='k',
                       node_size=250,
                       alpha=0.7
                       )
    nx.draw_networkx_nodes(G,pos,
                       nodelist=c2,
                       node_color='b',
                       node_size=250,
                       alpha=0.7
                       )
    nx.draw_networkx_nodes(G,pos,
                       nodelist=c1,
                       node_color='g',
                       node_size=250,
                       alpha=0.7
                       )
    plt.title('Raw Erdos-Reyni Graph with color-coded overlay of kmeans clustering k=3')
    nx.draw(G)
    plt.show()
    
    
