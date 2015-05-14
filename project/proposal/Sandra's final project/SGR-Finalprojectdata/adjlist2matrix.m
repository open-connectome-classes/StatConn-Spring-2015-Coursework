function [B]=adjlist2matrix(A)
%Returns the matrix of any adjacency list where the 
%input is A=any adj list with 2 columns 

rows = A(:,1);
cols = A(:,2);
s=ones(size(A,1),1);
BA = sparse([rows; cols],[cols; rows],[s,s]);
B=full(BA);
end