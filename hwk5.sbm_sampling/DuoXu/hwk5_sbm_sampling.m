%{

You should sample a graph from some model with some size, number of 
vertices, and some clusters

compute likelihood for same model with true number of clusters
compute likelihood for same model with some other number of clusters
report on the results (pictures, words, code)

%}


% sample
n=100;
p=0.21;
A=rand(n)<p;

% log-likelihood for SBM with K=1
loglikA=A*log(p)+(1-A)*log(1-p);
sumloglikA=sum(loglikA(:));


%% Creating Model and Sampling from It

numV = 100; % #vertices

sbmK_true = 3; % #classes/clusters
sbmB_true = SymmetricRandMatrix(sbmK_true); % Generate B of SBM
labelV_true = sort(floor(rand(numV,1) * sbmK_true) + 1); % uniformly labeling vertices with #sbmK classes

weightE_true = zeros(numV); % Construct distributions of each edges
for i = 1:sbmK_true
    for j = 1:sbmK_true
        weightE_true(find(labelV_true == i),find(labelV_true == j)) = sbmB_true(i,j);
    end
end

figure(1);
imagesc(weightE_true);
colormap gray, axis equal tight;

adjE = SymmetricRandMatrix(numV);
adjE = adjE < weightE_true;

figure(2);
imagesc(adjE);
colormap gray, axis equal tight;


%% Log Likelihood with Same(True) Cluster Number

labelV_same = kmeans(adjE, sbmK_true);
sbmB_same = zeros(sbmK_true);
for i = 1:sbmK_true
    for j = 1:sbmK_true
        sbmB_same(i,j) = mean(mean(adjE(find(labelV_same == i),find(labelV_same == j))));
    end
end

weightE_same = zeros(numV); % Construct distributions of each edges
for i = 1:sbmK_true
    for j = 1:sbmK_true
        weightE_same(find(labelV_same == i),find(labelV_same == j)) = sbmB_same(i,j);
    end
end

figure(3);
imagesc(weightE_same);
colormap gray, axis equal tight;

loglik_same = sum(sum(adjE * log(weightE_same) + (1-adjE) * log(1-weightE_same)));


%% Log Likelihood with Different Cluster Number

sbmK_diff = sbmK_true + 1;

labelV_diff = kmeans(adjE, sbmK_diff);
sbmB_diff = zeros(sbmK_diff);
for i = 1:sbmK_diff
    for j = 1:sbmK_diff
        sbmB_diff(i,j) = mean(mean(adjE(find(labelV_diff == i),find(labelV_diff == j))));
    end
end

weightE_diff = zeros(numV); % Construct distributions of each edges
for i = 1:sbmK_diff
    for j = 1:sbmK_diff
        weightE_diff(find(labelV_diff == i),find(labelV_diff == j)) = sbmB_diff(i,j);
    end
end

figure(4);
imagesc(weightE_diff);
colormap gray, axis equal tight;

loglik_diff = sum(sum(adjE * log(weightE_diff) + (1-adjE) * log(1-weightE_diff)));




