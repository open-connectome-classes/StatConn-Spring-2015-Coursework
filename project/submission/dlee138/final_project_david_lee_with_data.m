%
clear all;
%import data
url='https://www.dropbox.com/sh/idt3d0gylplyo31/AACqNPXHHbxKfYXuYXJ7un96a/dlee138new.zip?dl=0';
cmd=['wget ' url ' -O "dlee_138new.zip" --no-check-certificate'];
system(cmd);
unzip('dlee_138new.zip');
addpath('dlee_138new');

graph=importdata('Time1.mat');
graph2=importdata('Time2.mat');
graph3=importdata('Time3.mat');

figure(1);
imagesc(graph);
colormap(gray);
title('Figure 1: Graph at Time=1');


figure(2);
imagesc(graph2);
colormap(gray);
title('Figure 2: Graph at Time=2');


figure(3);
imagesc(graph3);
colormap(gray);
title('Figure 3: Graph at Time=3');


% Algorithm with 1 cluster
[likelihood_1,val_1] = fmincon(@(p) -sum(sum(graph*log(p)+(1-graph)*log(1-p))),0.01,[],[],[],[],0,1);
block1_t1=val_1;

% Algorithm with 2 cluster
block2_t1=1:10;
for i=1:10;
graph_2 = kmeans(graph,2);
clust2_11 = graph(graph_2==1,graph_2==1);
clust2_12 = graph(graph_2==1,graph_2==2);
clust2_22 = graph(graph_2==2,graph_2==2);

[likelihood_2,val_2] = fmincon(@(p) -sum(sum(clust2_11*log(p(1))+(1-clust2_11)*log(1-p(1))))-sum(sum(clust2_22*log(p(2))+(1-clust2_22)*log(1-p(2))))-sum(sum(clust2_12*log(p(3))+(1-clust2_12)*log(1-p(3)))),[0.01 0.01 0.001],[],[],[],[],[0 0 0],[1 1 1]);

block2_t1(i)=val_2;
end
val_2=mean(block2_t1);
%  Algorithm with 3 clusters (which is the "true" number of clusters in
%  this graph)
block3_t1=1:10;
for i=1:10;
graph_3 = kmeans(graph,3);

clust3_11 = graph(graph_3==1,graph_3==1);
clust3_12 = graph(graph_3==1,graph_3==2);
clust3_13 = graph(graph_3==1,graph_3==3);
clust3_22 = graph(graph_3==2,graph_3==2);
clust3_23 = graph(graph_3==2,graph_3==3);
clust3_33 = graph(graph_3==3,graph_3==3);

[likelihood_3,val_3] = fmincon(@(p) -sum(sum(clust3_11*log(p(1))+(1-clust3_11)*log(1-p(1))))-...
    sum(sum(clust3_22*log(p(2))+(1-clust3_22)*log(1-p(2))))-sum(sum(clust3_33*log(p(3))+(1-clust3_33)*log(1-p(3))))-sum(sum(clust3_12*log(p(4))+(1-clust3_12)*log(1-p(4))))-sum(sum(clust3_13*log(p(5))+(1-clust3_13)*log(1-p(5))))-sum(sum(clust3_23*log(p(6))+(1-clust3_23)*log(1-p(6)))),[0.01 0.01 0.01 0.001 0.001 0.001],[],[],[],[],[0 0 0 0 0 0],[1 1 1 1 1 1]);
block3_t1(i)=val_3;
end
val_3=mean(block3_t1);
% Figure 4
figure(4)
plot(1:3,-[val_1 val_2 val_3]);
xlabel('# of Clusters');
ylabel('Log Likelihood');
title('Figure 4: Log Likelihood vs. # of Clusters at Time=1');

% Algorithm with 1 cluster
[likelihood_1,val_1] = fmincon(@(p) -sum(sum(graph2*log(p)+(1-graph2)*log(1-p))),0.01,[],[],[],[],0,1)
block1_t2=val_1;

% Algorithm with 2 cluster
block2_t2=1:10;
for i=1:10;
graph2_2 = kmeans(graph2,2);
clust2_11 = graph2(graph2_2==1,graph2_2==1);
clust2_12 = graph2(graph2_2==1,graph2_2==2);
clust2_22 = graph2(graph2_2==2,graph2_2==2);

[likelihood_2,val_2] = fmincon(@(p) -sum(sum(clust2_11*log(p(1))+(1-clust2_11)*log(1-p(1))))-sum(sum(clust2_22*log(p(2))+(1-clust2_22)*log(1-p(2))))-sum(sum(clust2_12*log(p(3))+(1-clust2_12)*log(1-p(3)))),[0.01 0.01 0.001],[],[],[],[],[0 0 0],[1 1 1]);
block2_t2(i)=val_2;
end
val_2=mean(block2_t2);
%  Algorithm with 3 clusters (which is the "true" number of clusters in
%  this graph2)
block3_t2=1:10;
for i=1:10;
graph2_3 = kmeans(graph2,3);

clust3_11 = graph2(graph2_3==1,graph2_3==1);
clust3_12 = graph2(graph2_3==1,graph2_3==2);
clust3_13 = graph2(graph2_3==1,graph2_3==3);
clust3_22 = graph2(graph2_3==2,graph2_3==2);
clust3_23 = graph2(graph2_3==2,graph2_3==3);
clust3_33 = graph2(graph2_3==3,graph2_3==3);

[likelihood_3,val_3] = fmincon(@(p) -sum(sum(clust3_11*log(p(1))+(1-clust3_11)*log(1-p(1))))-...
    sum(sum(clust3_22*log(p(2))+(1-clust3_22)*log(1-p(2))))-sum(sum(clust3_33*log(p(3))+(1-clust3_33)*log(1-p(3))))-sum(sum(clust3_12*log(p(4))+(1-clust3_12)*log(1-p(4))))-sum(sum(clust3_13*log(p(5))+(1-clust3_13)*log(1-p(5))))-sum(sum(clust3_23*log(p(6))+(1-clust3_23)*log(1-p(6)))),[0.01 0.01 0.01 0.001 0.001 0.001],[],[],[],[],[0 0 0 0 0 0],[1 1 1 1 1 1]);
 block3_t2(i)=val_3;
end
val_3=mean(block3_t2);
% Figure 5
figure(5)
plot(1:3,-[val_1 val_2 val_3]);
xlabel('# of Clusters');
ylabel('Log Likelihood');
title('Figure 2: Log Likelihood vs. # of Clusters at Time=2');

% Algorithm with 1 cluster
block1_t3=1:10;
for i=1:10;
[likelihood_1,val_1] = fmincon(@(p) -sum(sum(graph3*log(p)+(1-graph3)*log(1-p))),0.01,[],[],[],[],0,1);
block1_t3(i)=val_1;
end
val_1=mean(block1_t3);
% Algorithm with 2 cluster
block2_t3=1:10;
for i=1:10;
graph3_2 = kmeans(graph3,2);
clust2_11 = graph3(graph3_2==1,graph3_2==1);
clust2_12 = graph3(graph3_2==1,graph3_2==2);
clust2_22 = graph3(graph3_2==2,graph3_2==2);

[likelihood_2,val_2] = fmincon(@(p) -sum(sum(clust2_11*log(p(1))+(1-clust2_11)*log(1-p(1))))-sum(sum(clust2_22*log(p(2))+(1-clust2_22)*log(1-p(2))))-sum(sum(clust2_12*log(p(3))+(1-clust2_12)*log(1-p(3)))),[0.01 0.01 0.001],[],[],[],[],[0 0 0],[1 1 1]);
block2_t3(i)=val_2;
end
val_2=mean(block2_t3);
%  Algorithm with 3 clusters (which is the "true" number of clusters in
%  this graph3)

block3_t3=1:100;
for i=1:100;
graph3_3 = kmeans(graph3,3);

clust3_11 = graph3(graph3_3==1,graph3_3==1);
clust3_12 = graph3(graph3_3==1,graph3_3==2);
clust3_13 = graph3(graph3_3==1,graph3_3==3);
clust3_22 = graph3(graph3_3==2,graph3_3==2);
clust3_23 = graph3(graph3_3==2,graph3_3==3);
clust3_33 = graph3(graph3_3==3,graph3_3==3);

[likelihood_3,val_3] = fmincon(@(p) -sum(sum(clust3_11*log(p(1))+(1-clust3_11)*log(1-p(1))))-...
    sum(sum(clust3_22*log(p(2))+(1-clust3_22)*log(1-p(2))))-sum(sum(clust3_33*log(p(3))+(1-clust3_33)*log(1-p(3))))-sum(sum(clust3_12*log(p(4))+(1-clust3_12)*log(1-p(4))))-sum(sum(clust3_13*log(p(5))+(1-clust3_13)*log(1-p(5))))-sum(sum(clust3_23*log(p(6))+(1-clust3_23)*log(1-p(6)))),[0.01 0.01 0.01 0.001 0.001 0.001],[],[],[],[],[0 0 0 0 0 0],[1 1 1 1 1 1]);
block3_t3(i)=val_3;
end
val_3=mean(block3_t3);
 
% Figure 3
figure(6)
plot(1:3,-[val_1 val_2 val_3]);
xlabel('# of Clusters');
ylabel('Log Likelihood');
title('Figure 6: Log Likelihood vs. # of Clusters at Time=6');

predicted_blocks_t1=1;
lik_t1=[-mean(block1_t1),-mean(block2_t1),-mean(block3_t1)]
for i=1:2
    if (lik_t1(i)<1.2*lik_t1(i+1))
        predicted_blocks_t1=predicted_blocks_t1+1;
    end
end
predicted_blocks_t1

predicted_blocks_t2=1;
lik_t2=[-mean(block1_t2),-mean(block2_t2),-mean(block3_t2)]
for i=1:2
    if (lik_t2(i)<1.2*lik_t2(i+1))
        predicted_blocks_t2=predicted_blocks_t2+1;
    end
end
predicted_blocks_t2

predicted_blocks_t3=1;
lik_t3=[-mean(block1_t3),-mean(block2_t3),-mean(block3_t3)]
for i=1:2
    if (lik_t3(i)<1.2*lik_t3(i+1))
        predicted_blocks_t3=predicted_blocks_t3+1;
    end
end
predicted_blocks_t3
        
figure(7)
predictions=[predicted_blocks_t1,predicted_blocks_t2,predicted_blocks_t3]
bar(predictions)
title('Predicted number of blocks at each time point');
ylabel('Number of Predicted Blocks');
xlabel('Time point');