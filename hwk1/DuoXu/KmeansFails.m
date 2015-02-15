

numNode = 6;
A = rand(numNode, numNode);
threshold = 0.7;
A(A > threshold) = 1;
A(A ~= 1) = 0;
A = triu(A);
coordinates = rand(numNode,2) * 10;

ind = kmeans(A', 2);
group1 = find(ind==1);
group2 = find(ind==2);

figure
hold on
gplot(A, coordinates, '-*');
scatter(coordinates(group1,1), coordinates(group1,2), 'r');
hold off

