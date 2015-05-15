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
