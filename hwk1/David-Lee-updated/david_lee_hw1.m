
%David Lee
%Statistical Connectomics HW1

%adjacency matrix 
a = [0,1,1,0,0,0;1,0,1,0,0,0;1,1,0,1,0,0;0,0,1,0,1,1;0,0,0,1,0,1;0,0,0,1,1,0];
clusters = kmeans(a, 2)

%Plots out the graph which consists of 6 vertices and 7 egdes
Coordinates = [0,1;1,0;1,1;2,2;2,3;3,2];
gplot(a,Coordinates)
hold all

%Visually represents results of k-means clustering by seperating verticies into two categories (Red and Green)
for i=1:6
if clusters(i,1)==1
scatter(Coordinates(i,1),Coordinates(i,2),'g','filled')
end
end
for i=1:6
if clusters(i,1)==2
scatter(Coordinates(i,1),Coordinates(i,2),'r','filled')
end
end