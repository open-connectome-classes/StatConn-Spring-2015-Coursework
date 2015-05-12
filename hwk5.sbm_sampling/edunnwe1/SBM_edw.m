% You should sample a graph from some model with some size, number of 
% vertices, and some clusters
% compute likelihood for same model with true number of clusters
% compute likelihood for same model with some other number of clusters
% report on the results (pictures, words, code)


%% SBM

n = 100; % #vertices

K1 = 3; % true number of classes/clusters
B = triu(rand(K1));
B = B + tril(B',-1); %symmetric block matrix with probabilites for an edge within (on the diagonal) and between (off diagonal) clusters

clusters = sort(floor(rand(n,1) * K1) + 1); % cluster vertices

probability_edge = zeros(n); % sample a graph given known clusters of vertices
for i = 1:K1
    for j = 1:K1
        probability_edge(clusters == i,clusters == j) = B(i,j);
    end
end

A = triu(rand(n));
A = A + tril(A',-1);
A = A < probability_edge; %sample the adjacency matrix according to B

figure('Position',[500 300 800 300]);
subplot(1,2,1);
imagesc(probability_edge);
colormap gray, axis equal tight; 
title('B');colorbar


subplot(1,2,2);
imagesc(A);
colormap gray, axis equal tight;
title('Adjacency matrix, sorted into clustered blocks');colorbar


%% Log LiK1elihood with with true number of clusters

loglikA_K1=A*log(probability_edge)+(1-A)*log(1-probability_edge);
sumloglikA_K1=sum(loglikA(:))

%% Log LiK1elihood with false number of clusters

K0 = K1 + 1; %false number of clusters
%(notation: 0 for false, 1 for true);

B0 = triu(rand(K0));
B0 = B0 + tril(B0',-1); %symmetric block matrix with probabilites for an edge within (on the diagonal) and between (off diagonal) clusters

clusters0 = sort(floor(rand(n,1) * K0) + 1); % cluster vertices

probability_edge0 = zeros(n); % Construct distributions of each edges
for i = 1:K0
    for j = 1:K0
        probability_edge0(clusters0 == i,clusters0 == j) = B0(i,j);
    end
end



figure('Position',[200 300 1200 300]);
subplot(1,3,1);
imagesc(probability_edge);
colormap gray, axis equal tight; 
title('True B');colorbar

subplot(1,3,2);
imagesc(probability_edge0);
colormap gray, axis equal tight; 
title('B based on false number of clusters');colorbar

subplot(1,3,3);
imagesc(A);
colormap gray, axis equal tight;
title('Adjacency matrix, sorted into true clustered blocks');colorbar
