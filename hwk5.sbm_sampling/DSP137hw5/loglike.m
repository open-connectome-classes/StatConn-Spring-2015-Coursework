function [sloglikA] = loglike(A, maxk)
% Input:  A    = adjacency matrix for undirected graph
%         maxk = maximum number of clusters
% Output: sloglikA = vector of sum-log-likelihood values for each
%                    clustering of A
% HGP 2015 
% Statistical Connectomics HW05
%
% S.D.G.

%maxk = 10;
%esumloglikA = NaN(maxk,1); 
%meanloglikA = NaN(maxk,1);
sloglikA = NaN(maxk,1); % initializes sloglikA
x = [1:maxk]';  
n=length(A);
% This will create a 'mock' probability matrix based on the proportions of
% edges in A based on the clustering given by kmeans for 'c' clusters.
% Following the creation of 'Bhat' we calculate the log-likelihood for each
% entry in A. 
for c = 1:maxk
    clust = kmeans(A,c);
    Bp = NaN(c,c);
    for i=1:c
        for j=1:c
            Bp(i,j)=sum(sum(A(clust==i,clust==j)))/(sum(clust==i)*sum(clust==j));
        end
    end
    likA = NaN(length(A));
    loglikA = NaN(length(A));
    for i=1:n
        for j=1:n
            p = Bp(clust(i),clust(j));
            a = A(i,j);
            likA(i,j)=p^a*(1-p)^(1-a);
            loglikA(i,j)=a*log(p)+(1-a)*log(1-p);
        end
    end
    prodlikA=prod(likA(:)); % produces the product of the loglikelihood
    sumloglikA=sum(loglikA(:));

    % esumloglikA(c) = exp(sumloglikA);
    % meanloglikA(c) = -mean(sumloglikA)/(n^2);
    sloglikA(c) = sumloglikA; %keeps track of sum-log-like(A) based on 
                              %number of clusters used
end

%Plots the sum-log-likelihood against number of clusters
h2 = figure;
plot(x,sloglikA,'black-')
xlabel('number of clusters')
ylabel('log-likelihood')
print(h2,'-djpeg','hw5_logl')
end
