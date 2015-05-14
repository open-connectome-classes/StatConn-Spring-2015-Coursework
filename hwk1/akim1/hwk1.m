m = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];
imagesc(m);
title('adjacency matrix (result of kmeans = [2 1 2 1])');
kmeans(m,2)
