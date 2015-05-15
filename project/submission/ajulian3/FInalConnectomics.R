# Final Project
# Statistical Connectomics
# Alan Juliano


install.packages('e1071',dependencies=TRUE) 
install.packages("MASS")
install.packages("brainwaver")
library(e1071)
library(MASS)
install.packages("class")
install.packages("DAAG")
library(class)
library("DAAG")
library("brainwaver")
data(brain)
data<-as.matrix(brain)


# In our data, V8-V15 and V22-30 are the regions of interest when examining ADHD
superiorpc<-data[,8:15]
premotor<-data[,22:30]

#Depression effects activity in the V56-V60 (Hippocampus), V63-65(Amygdala), V32-35(Striatum), V78-82 #(Thalamus)  

hippo<-data[,56:60]
amyg<-data[,63:65]
striatum<-data[,32:35]
thalamus<-data[,78:82]

# Correlation Matrices for each level of wavelet decomposition
#ADHD Regions

wave.cor.list1<-const.cor.list(superiorpc, method = "modwt" ,wf = "la8", n.levels = 6,
boundary = "periodic", p.corr = 0.975)

wave.cor.list2<-const.cor.list(premotor, method = "modwt" ,wf = "la8", n.levels = 6,
 boundary = "periodic", p.corr = 0.975)
  
#Depression Regions

wave.cor.list3<-const.cor.list(hippo, method = "modwt" ,wf = "la8", n.levels = 6,
boundary = "periodic", p.corr = 0.975)

wave.cor.list4<-const.cor.list(amyg, method = "modwt" ,wf = "la8", n.levels = 6,
boundary = "periodic", p.corr = 0.975)

wave.cor.list5<-const.cor.list(striatum, method = "modwt" ,wf = "la8", n.levels = 6,
boundary = "periodic", p.corr = 0.975)

wave.cor.list6<-const.cor.list(thalamus, method = "modwt" ,wf = "la8", n.levels = 6,
boundary = "periodic", p.corr = 0.975)

# Adjacency Matrix within each of our partitioned regions along with required edges for a graph to be generated
#nbedges=how many edges are required for a graph to be generated
#sup.thresh= threshold for each correlation matrix to generate an adgacency matrix with correctly defined edges

# ADHD Region of Brain
#Region 1
adj.mat.list1<-const.adj.mat(wave.cor.list1[[4]], sup = 0.44,proc.length=dim(data)[1],
num.levels=6)
nb.edges1<-sum(adj.mat.list1)/2
sup.thresh1<-choose.thresh.nbedges(wave.cor.list1[[4]],nb.edges=nb.edges1, proc.length=dim(data)[1],num.levels=6)
pvalue.cor1<-p.value.compute(wave.cor.list1[[4]],proc.length=dim(data)[1], sup=0.44, num.levels=6)
pvalue.thresh1<-compute.FDR(pvalue.cor1,.05)


#Region 2
adj.mat.list2<-const.adj.mat(wave.cor.list2[[4]], sup = 0.44,proc.length=dim(data)[1],
num.levels=6)
nb.edges2<-sum(adj.mat.list2)/2
sup.thresh2<-choose.thresh.nbedges(wave.cor.list2[[4]],nb.edges=nb.edges2, proc.length=dim(data)[1],num.levels=6)
pvalue.cor2<-p.value.compute(wave.cor.list2[[4]],proc.length=dim(data)[1], sup=0.44, num.levels=6)
pvalue.thresh2<-compute.FDR(pvalue.cor2,.05)

#Depression Region of Brain
#Region 3
adj.mat.list3<-const.adj.mat(wave.cor.list3[[4]], sup = 0.44,proc.length=dim(data)[1],
num.levels=6)
nb.edges3<-sum(adj.mat.list3)/2
sup.thresh3<-choose.thresh.nbedges(wave.cor.list3[[4]],nb.edges=nb.edges3, proc.length=dim(data)[1],num.levels=6)
pvalue.cor3<-p.value.compute(wave.cor.list3[[4]],proc.length=dim(data)[1], sup=0.44, num.levels=6)
pvalue.thresh3<-compute.FDR(pvalue.cor3,.05)
#Region 4
adj.mat.list4<-const.adj.mat(wave.cor.list4[[4]], sup = 0.44,proc.length=dim(data)[1],
num.levels=6)
nb.edges4<-sum(adj.mat.list4)/2
sup.thresh4<-choose.thresh.nbedges(wave.cor.list4[[4]],nb.edges=nb.edges4, proc.length=dim(data)[1],num.levels=6)
pvalue.cor4<-p.value.compute(wave.cor.list4[[4]],proc.length=dim(data)[1], sup=0.44, num.levels=6)
pvalue.thresh4<-compute.FDR(pvalue.cor4,.05)
#Region5
adj.mat.list5<-const.adj.mat(wave.cor.list5[[4]], sup = 0.44,proc.length=dim(data)[1],
num.levels=6)
nb.edges5<-sum(adj.mat.list5)/2
sup.thresh5<-choose.thresh.nbedges(wave.cor.list5[[4]],nb.edges=nb.edges5, proc.length=dim(data)[1],num.levels=6)
pvalue.cor5<-p.value.compute(wave.cor.list5[[4]],proc.length=dim(data)[1], sup=0.44, num.levels=6)
pvalue.thresh5<-compute.FDR(pvalue.cor5,.05)
#Region 6
adj.mat.list6<-const.adj.mat(wave.cor.list6[[4]], sup = 0.44,proc.length=dim(data)[1],
num.levels=6)
nb.edges6<-sum(adj.mat.list6)/2
sup.thresh6<-choose.thresh.nbedges(wave.cor.list6[[4]],nb.edges=nb.edges6, proc.length=dim(data)[1],num.levels=6)

# False Discovery Rate
# designed to control the expected proportion of incorrectly rejected null hypotheses (or else large type 1 errors)
pvalue.thresh1<-compute.FDR(pvalue.cor1,.05)
pvalue.thresh2<-compute.FDR(pvalue.cor2,.05)
pvalue.thresh3<-compute.FDR(pvalue.cor3,.05)
pvalue.thresh4<-compute.FDR(pvalue.cor4,.05)
pvalue.thresh5<-compute.FDR(pvalue.cor5,.05)
pvalue.thresh6<-compute.FDR(pvalue.cor6,.05)
# pvalue.thresh1 and pvalue.thresh3 have the highest p values for FDR 

#Correlations Between Adjacency Matrix and Correlations
adj.mat1<-correlations.to.adjacencies(wave.cor.list1,edge.func=(function(x){x*log(x)}))
adj.mat2<-correlations.to.adjacencies(wave.cor.list2,edge.func=(function(x){x*log(x)}))
adj.mat3<-correlations.to.adjacencies(wave.cor.list3,edge.func=(function(x){x*log(x)}))
adj.mat4<-correlations.to.adjacencies(wave.cor.list4,edge.func=(function(x){x*log(x)}))
adj.mat5<-correlations.to.adjacencies(wave.cor.list5,edge.func=(function(x){x*log(x)}))
adj.mat6<-correlations.to.adjacencies(wave.cor.list6,edge.func=(function(x){x*log(x)}))

#Determining Graph Efficiency and Compare between ADHD and Depression Regions



n.regions<-dim(brain)[2]

#Construction of the correlation matrices for each level of the wavelet decomposition
wave.cor.list<-const.cor.list(brain, method = "modwt" ,wf = "la8", n.levels = 6, 
                               boundary = "periodic", p.corr = 0.975)
# First we must implement a threshold on correlation
supseq1<-((1:10)/10) #sequence of the correlation threshold 
gmax<-length(supseq1)
Eglob<-matrix(0,6,gmax)
Eloc<-matrix(0,6,gmax)
Cost<-matrix(0,6,gmax)

n.levels<-6

# A for loop was then generated 

for(i in 1:gmax){
n.sup<-supseq1[i]

#Construction of the adjacency matrices associated to each level of the wavelet decomposition
wave.adj.list<-const.adj.list(wave.cor.list, sup = n.sup)

#For each level of the wavelet decomposition
for(j in 1:n.levels){

Eglob.brain<-global.efficiency(wave.adj.list[[j]],
                        weight.mat=matrix(1,n.regions,n.regions))
Eglob[j,i]<-Eglob.brain$eff

Eloc.brain<-local.efficiency(wave.adj.list[[j]],
                                weight.mat=matrix(1,n.regions,n.regions))
Eloc[j,i]<-Eloc.brain$eff

Cost.brain<-global.cost(wave.adj.list[[j]],
                                weight.mat=matrix(1,n.regions,n.regions))
Cost[j,i]<-Cost.brain

}}

plot(sup.seq,(1:gmax)/2,type='n',xlab='R: Correlation Threshold',ylab='Global Efficiency (Graph)',
     cex.axis=2,cex.lab=2,xlim=c(0,1),ylim=c(0,1))

for(i in 1:n.levels){
lines(sup.seq,Eglob[i,],type='l',col=i,lwd=2)
}
legend(x="topright",legend=c("Level 1","Level 2","Level 3","Level 4",
                                "Level 5","Level 6"),fill=TRUE,col=(1:n.levels),lwd=2)

