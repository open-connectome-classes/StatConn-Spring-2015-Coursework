%
clear all;
c1 = 30;
c2 = 30;
c3 = 30;
p_11 = 0.7;
p_12 = 0.05;
p_13 = 0.05;
p_21 = 0.05;
p_22 = 0.7;
p_23 = 0.05;
p_31 = 0.05;
p_32 = 0.05;
p_33 = 0.7;

graph = zeros(c1+c2+c3);

% SBM
graph(1:c1,1:c1) = rand(c1)<p_11;
graph(1:c1,c1+1:c1+c2) = rand(c1,c2)<p_12;
graph(1:c1,c1+c2+1:c1+c2+c3) = rand(c1,c3)<p_13;
X = triu(graph(1:c1,1:c1))+transpose(triu(graph(1:c1,1:c1),1));
graph(1:c1,1:c1)=X;

graph(c1+1:c1+c2,1:c1) = transpose(graph(1:c1,c1+1:c1+c2)); 
graph(c1+1:c1+c2,c1+1:c1+c2) =  rand(c2)<p_22;
graph(c1+1:c1+c2,c1+c2+1:c1+c2+c3) = rand(c2,c3)<p_23;
Y = triu(graph(c1+1:c1+c2,c1+1:c1+c2))+transpose(triu(graph(c1+1:c1+c2,c1+1:c1+c2),1));
graph(c1+1:c1+c2,c1+1:c1+c2)=Y;

graph(c1+c2+1:c1+c2+c3,1:c1) = transpose(graph(1:c1,c1+c2+1:c1+c2+c3));
graph(c1+c2+1:c1+c2+c3,c1+1:c1+c2) = transpose(graph(c1+1:c1+c2,c1+c2+1:c1+c2+c3));
graph(c1+c2+1:c1+c2+c3,c1+c2+1:c1+c2+c3) = rand(c3)<p_33;
Z = triu(graph(c1+c2+1:c1+c2+c3,c1+c2+1:c1+c2+c3))+transpose(triu(graph(c1+c2+1:c1+c2+c3,c1+c2+1:c1+c2+c3),1));
graph(c1+c2+1:c1+c2+c3,c1+c2+1:c1+c2+c3)=Z;

graph2 = graph;
p_11 = 0.7;
p_12 = 0.7;
p_21 = 0.7;
p_22 = 0.7;

graph2(1:c1,1:c1) = rand(c1)<p_11;
graph2(1:c1,c1+1:c1+c2) = rand(c1,c2)<p_12;
X = triu(graph2(1:c1,1:c1))+transpose(triu(graph2(1:c1,1:c1),1));
graph2(1:c1,1:c1)=X;

graph2(c1+1:c1+c2,1:c1) = transpose(graph2(1:c1,c1+1:c1+c2));
graph2(c1+1:c1+c2,c1+1:c1+c2) = rand(c2)<p_22;
Y = triu(graph2(c1+1:c1+c2,c1+1:c1+c2))+transpose(triu(graph2(c1+1:c1+c2,c1+1:c1+c2),1));
graph2(c1+1:c1+c2,c1+1:c1+c2)=Y;

graph3 = graph2;
p_11 = 0.05;
p_12 = 0.7;
p_21 = 0.7;
p_22 = 0.05;

graph3(1:c1,1:c1) = rand(c1)<p_11;
graph3(1:c1,c1+1:c1+c2) = rand(c1,c2)<p_12;
X = triu(graph3(1:c1,1:c1))+transpose(triu(graph3(1:c1,1:c1),1));
graph3(1:c1,1:c1)=X;


graph3(c1+1:c1+c2,1:c1) = transpose(graph3(1:c1,c1+1:c1+c2));
graph3(c1+1:c1+c2,c1+1:c1+c2) = rand(c2)<p_22;
Y = triu(graph3(c1+1:c1+c2,c1+1:c1+c2))+transpose(triu(graph3(c1+1:c1+c2,c1+1:c1+c2),1));
graph3(c1+1:c1+c2,c1+1:c1+c2)=Y;

% Figure 1
figure(1)
imagesc(graph);
colormap(gray);
figure(2)
imagesc(graph2);
colormap(gray);
figure(3)
imagesc(graph3);
colormap(gray);


% Algorithm with 1 cluster
[likelihood_1,val_1] = fmincon(@(p) -sum(sum(graph*log(p)+(1-graph)*log(1-p))),0.01,[],[],[],[],0,1)

% Algorithm with 2 cluster
graph_2 = kmeans(graph,2);
clust2_11 = graph(graph_2==1,graph_2==1);
clust2_12 = graph(graph_2==1,graph_2==2);
clust2_22 = graph(graph_2==2,graph_2==2);

[likelihood_2,val_2] = fmincon(@(p) -sum(sum(clust2_11*log(p(1))+(1-clust2_11)*log(1-p(1))))-sum(sum(clust2_22*log(p(2))+(1-clust2_22)*log(1-p(2))))-sum(sum(clust2_12*log(p(3))+(1-clust2_12)*log(1-p(3)))),[0.01 0.01 0.001],[],[],[],[],[0 0 0],[1 1 1])

%  Algorithm with 3 clusters (which is the "true" number of clusters in
%  this graph)
graph_3 = kmeans(graph,3);

clust3_11 = graph(graph_3==1,graph_3==1);
clust3_12 = graph(graph_3==1,graph_3==2);
clust3_13 = graph(graph_3==1,graph_3==3);
clust3_22 = graph(graph_3==2,graph_3==2);
clust3_23 = graph(graph_3==2,graph_3==3);
clust3_33 = graph(graph_3==3,graph_3==3);

[likelihood_3,val_3] = fmincon(@(p) -sum(sum(clust3_11*log(p(1))+(1-clust3_11)*log(1-p(1))))-...
    sum(sum(clust3_22*log(p(2))+(1-clust3_22)*log(1-p(2))))-sum(sum(clust3_33*log(p(3))+(1-clust3_33)*log(1-p(3))))-sum(sum(clust3_12*log(p(4))+(1-clust3_12)*log(1-p(4))))-sum(sum(clust3_13*log(p(5))+(1-clust3_13)*log(1-p(5))))-sum(sum(clust3_23*log(p(6))+(1-clust3_23)*log(1-p(6)))),[0.01 0.01 0.01 0.001 0.001 0.001],[],[],[],[],[0 0 0 0 0 0],[1 1 1 1 1 1])

% Figure 2
figure(4)
plot(1:3,-[val_1 val_2 val_3]);
xlabel('# of Clusters');
ylabel('Log Likelihood');
title('Figure 2: Log Likelihood vs. # of Clusters');

% Algorithm with 1 cluster
[likelihood_1,val_1] = fmincon(@(p) -sum(sum(graph2*log(p)+(1-graph2)*log(1-p))),0.01,[],[],[],[],0,1)

% Algorithm with 2 cluster
graph2_2 = kmeans(graph2,2);
clust2_11 = graph2(graph2_2==1,graph2_2==1);
clust2_12 = graph2(graph2_2==1,graph2_2==2);
clust2_22 = graph2(graph2_2==2,graph2_2==2);

[likelihood_2,val_2] = fmincon(@(p) -sum(sum(clust2_11*log(p(1))+(1-clust2_11)*log(1-p(1))))-sum(sum(clust2_22*log(p(2))+(1-clust2_22)*log(1-p(2))))-sum(sum(clust2_12*log(p(3))+(1-clust2_12)*log(1-p(3)))),[0.01 0.01 0.001],[],[],[],[],[0 0 0],[1 1 1])

%  Algorithm with 3 clusters (which is the "true" number of clusters in
%  this graph2)
graph2_3 = kmeans(graph2,3);

clust3_11 = graph2(graph2_3==1,graph2_3==1);
clust3_12 = graph2(graph2_3==1,graph2_3==2);
clust3_13 = graph2(graph2_3==1,graph2_3==3);
clust3_22 = graph2(graph2_3==2,graph2_3==2);
clust3_23 = graph2(graph2_3==2,graph2_3==3);
clust3_33 = graph2(graph2_3==3,graph2_3==3);

[likelihood_3,val_3] = fmincon(@(p) -sum(sum(clust3_11*log(p(1))+(1-clust3_11)*log(1-p(1))))-...
    sum(sum(clust3_22*log(p(2))+(1-clust3_22)*log(1-p(2))))-sum(sum(clust3_33*log(p(3))+(1-clust3_33)*log(1-p(3))))-sum(sum(clust3_12*log(p(4))+(1-clust3_12)*log(1-p(4))))-sum(sum(clust3_13*log(p(5))+(1-clust3_13)*log(1-p(5))))-sum(sum(clust3_23*log(p(6))+(1-clust3_23)*log(1-p(6)))),[0.01 0.01 0.01 0.001 0.001 0.001],[],[],[],[],[0 0 0 0 0 0],[1 1 1 1 1 1])
 
% Figure 2
figure(5)
plot(1:3,-[val_1 val_2 val_3]);
xlabel('# of Clusters');
ylabel('Log Likelihood');
title('Figure 2: Log Likelihood vs. # of Clusters');

% Algorithm with 1 cluster
[likelihood_1,val_1] = fmincon(@(p) -sum(sum(graph3*log(p)+(1-graph3)*log(1-p))),0.01,[],[],[],[],0,1)

% Algorithm with 2 cluster
graph3_2 = kmeans(graph3,2);
clust2_11 = graph3(graph3_2==1,graph3_2==1);
clust2_12 = graph3(graph3_2==1,graph3_2==2);
clust2_22 = graph3(graph3_2==2,graph3_2==2);

[likelihood_2,val_2] = fmincon(@(p) -sum(sum(clust2_11*log(p(1))+(1-clust2_11)*log(1-p(1))))-sum(sum(clust2_22*log(p(2))+(1-clust2_22)*log(1-p(2))))-sum(sum(clust2_12*log(p(3))+(1-clust2_12)*log(1-p(3)))),[0.01 0.01 0.001],[],[],[],[],[0 0 0],[1 1 1])

%  Algorithm with 3 clusters (which is the "true" number of clusters in
%  this graph3)
graph3_3 = kmeans(graph3,3);

clust3_11 = graph3(graph3_3==1,graph3_3==1);
clust3_12 = graph3(graph3_3==1,graph3_3==2);
clust3_13 = graph3(graph3_3==1,graph3_3==3);
clust3_22 = graph3(graph3_3==2,graph3_3==2);
clust3_23 = graph3(graph3_3==2,graph3_3==3);
clust3_33 = graph3(graph3_3==3,graph3_3==3);

[likelihood_3,val_3] = fmincon(@(p) -sum(sum(clust3_11*log(p(1))+(1-clust3_11)*log(1-p(1))))-...
    sum(sum(clust3_22*log(p(2))+(1-clust3_22)*log(1-p(2))))-sum(sum(clust3_33*log(p(3))+(1-clust3_33)*log(1-p(3))))-sum(sum(clust3_12*log(p(4))+(1-clust3_12)*log(1-p(4))))-sum(sum(clust3_13*log(p(5))+(1-clust3_13)*log(1-p(5))))-sum(sum(clust3_23*log(p(6))+(1-clust3_23)*log(1-p(6)))),[0.01 0.01 0.01 0.001 0.001 0.001],[],[],[],[],[0 0 0 0 0 0],[1 1 1 1 1 1])

 
% Figure 3
figure(6)
plot(1:3,-[val_1 val_2 val_3]);
xlabel('# of Clusters');
ylabel('Log Likelihood');
title('Figure 3: Log Likelihood vs. # of Clusters');

 