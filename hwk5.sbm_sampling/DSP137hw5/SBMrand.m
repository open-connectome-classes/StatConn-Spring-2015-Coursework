function [B,v,A] = SBMrand(n,k,s)
% Input:  n = number of vertices
%         k = number of clusters
%         s = if 0, don't sort vertices, else sort. 
%                default is to sort the vertices by cluster
% Output: B = block probability matrix
%         v = vector denoting to which cluster 
%             each vertex belongs
%         A = random matrix with k clusters
%         vp= vector denoting clusters to vertices sorted
% HGP 2015 
% Statistical Connectomics HW05
%
% S.D.G.

if (nargin<3)
    s=1;
end
    
B = NaN(k);
% %This will create a completely 
% %randomized symmetric block matrix
% for i=1:k
%     for j=1:i
%         B(i,j)=rand(1);
%         if (j~=i)
%             B(j,i)=B(i,j);
%         end
%     end
% end 

%This will create a block matrix depicting higher probabilities of an edge
%between vertices in the same cluster and lower probabilities of an edge
%between vertices in different clusters. 
for i=1:k
    for j=1:i
        if j==i
            B(i,j)=.75+(.25)*rand; % prob on diag are >= .75
        else 
            B(i,j)=(.25)*rand; % creates probs for off diag <= .25
            B(j,i)=B(i,j);  % forces symmetry
        end
    end
end 

% creates a vector of numbers between 1 and k as cluster assignments for
% vertices
v0 = NaN(n,1);
for i=1:n
    if i<=k
        v0(i)=i; % ensures at least 1 element per cluster
    else
        v0(i)=randi(k);
    end
end

% if not required to sort, vertices in clusters are spread out in matrix.
if s == 0 
    v=v0;
else
    v=sort(v0); % relabels vertices according to cluster
end

% To create random matrix A with probability distribution B and vertex
% assignment v, for each entry on upper diagonal generate a random number
% between 0 and 1. If said number is less than the probability given by the
% corresponding entry in B (based on vertex assignment), gives a 1 (edge)
% and otherwise a 0 (non-edge). Lower diagonal is generated as you go
% forcing symmetry (undirected graph).
A = NaN(n);
for i=1:n
    for j=1:i
        A(i,j)=rand<B(v(i),v(j)); 
        if (j~=i)
            A(j,i)=A(i,j);
        end
    end
end

        
