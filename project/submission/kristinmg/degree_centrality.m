function degree = degree_centrality(A)

% Convert A to zeros and ones (1 if there is an edge, else 0)
B = logical(A);

% Calculate degree
degree = sum(B,2);

        