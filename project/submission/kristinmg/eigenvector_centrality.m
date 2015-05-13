function eig_centrality = eigenvector_centrality(A)

[V,D] = eig(A);
[maxeig,indexmax] = max(diag(D));  % The maximum eigenvalue and its index
eig_centrality = abs(V(:,indexmax)); % The leading eigenvector