function w_degree = weighted_degree(A)

w_degree = zeros(size(A,1),1);
for i = 1:size(A,1)
    w_degree(i) = sum(A(i,:));
end