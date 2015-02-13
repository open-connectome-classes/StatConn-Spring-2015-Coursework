#
# EN.580.694 Statistical Connectomics
# HW01 Sp2015
# 
# Heather G. Patsolic
# S.D.G.

library(igraph)

x<-sample(0:1,100,replace=T)
A<-matrix(x,10,10)
A<-t(A)*A;
g<-graph.adjacency(A,mode="undirected",weighted=NULL,diag=T,add.colnames=NULL,add.rownames=NA)

plot(g)

kcluster<-kmeans(A,centers=2)
kcluster3<-kmeans(A,3)

B <- matrix(rbind(c(1,0,0,0,1,0,1,0,0,0),c(0,0,0,0,1,1,0,0,0,1),c(0,0,1,1,0,1,0,0,0,0),
c(0,0,1,0,0,0,0,0,0,0),c(1,1,0,0,1,0,0,0,1,0),c(0,1,1,0,0,1,0,0,1,0),
c(1,0,0,0,0,0,0,0,1,1),c(0,0,0,0,0,0,0,1,0,0),c(0,0,0,0,1,1,1,0,0,0),
c(0,1,0,0,0,0,1,0,0,0)),10,10)

graphB<-graph.adjacency(A,mode="undirected",weighted=NULL,diag=T,add.colnames=NULL,add.rownames=NA)
plot(graphB)

kb<-kmeans(B,2)
kb$cluster
#1 1 2 2 2 1 2 1 2 2 2
kb3<-kmeans(B,3)
kb3$cluster
#3 3 1 1 2 1 2 1 3 1




#
##
### Soli Deo Gloria
