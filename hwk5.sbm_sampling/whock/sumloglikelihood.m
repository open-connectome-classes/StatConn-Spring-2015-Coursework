function [sumLL] = sumloglikelihood(c0,ctest)
%% Loglikelihood fx draft version 1
% Will Hockeimer JHU Neuroscience
baseps = [0.75,0.05];
[adj,betamat] = sbm(100,c0,baseps,assortive);
adj_test = kmeans(adj,ctest);
sumLL=sum(adj_test*log(betamat)+(1-adj_test)*log(1-betamat));
end