function [adj,betamat] = sbm(nverts,kblocks,baseps,assortive)
%% EN.580.694 Statistical Connectomics
% Homework 5 - Sample from an SBM and compute log-likelihood
% William Hockeimer, JHU Neuroscience
% Parameters    nverts -- number of vertices
%               kblocks -- number of blocks (assume symmetric graph)
%               baseps -- high and low prob for on/off diag base prob
%               assortive -- pass t/f for whether graph is assortive

%%
if nargin <4
    baseps = fliplr(baseps); %if disassortive, 
end

verts_per_block = floor(nverts/kblocks);
labels=[];
for i = 1:kblocks
    labels = [labels,i*ones(1,verts_per_block)];
end

%% Defin e the beta matrix of inter- and intra-block connection probs
betamat = zeros(kblocks,kblocks);

for ii = 1:kblocks
    for jj = 1:kblocks
        if ii == jj
            betamat(ii,jj) = baseps(1) + 0.15*rand();
        else
            betamat(ii,jj) = baseps(2) + 0.1*rand();
        end
        
    end
end

%% Define adjacency matrix using label vector and beta matrix

adj = zeros(nverts,nverts);
for r = 1:nverts
    for c = 1:nverts
        adj(r,c) = rand() < betamat(labels(r),labels(c));    
    end
end

%% Plot
figure(1);
spy(adj)
title('Adjacency Matrix','FontSize',18)
xlabel('Vertex','FontSize',16)
ylabel('Vertex','FontSize',16)



end