function cosine_similarity = cossim(A)

cosine_similarity = zeros(size(A));
for i = 1:size(A,1)
    for j = 1:size(A,2)
        cosine_similarity(i,j) = dot(A(i,:),A(j,:))/(norm(A(i,:))*norm(A(j,:)));
    end
end