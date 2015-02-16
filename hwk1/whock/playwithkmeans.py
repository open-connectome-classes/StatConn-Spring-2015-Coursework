import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.vq import *
import networkx as nx
from sklearn.cluster import KMeans
from datetime import date

def playwithkmeans(n=50,k=3,p=0.6,save=False):
    '''randomly generate a random Watts-Strogatz graph with 
    n - nodes
    k - connected to k neighbors
    p - rewiring from base NN ring with prob b. 
    Labeled graph is plotted according to kmeans clust=3 
    to see by eye how "well" it does 
    WH StatConn HW 1
    '''
    G = nx.watts_strogatz_graph(n,k,p)
    pos = nx.random_layout(G)
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
    nx.draw_networkx_nodes(G,pos,nodelist=c1,node_color='r',node_size=500,alpha=0.8)
    nx.draw_networkx_nodes(G,pos,nodelist=c2,node_color='g',node_size=500,alpha=0.8)
    nx.draw_networkx_nodes(G,pos,nodelist=c3,node_color='b',node_size=500,alpha=0.8)
    nx.draw_networkx_edges(G,pos)
    plt.title('Random Graph G with color-coded overlay of kmeans clustering k=3')
    if save:
        plt.savefig('C:\\Users\Will\\Pictures\\graph-%s.pdf'%date.today(),format='pdf')
    plt.show()
    
