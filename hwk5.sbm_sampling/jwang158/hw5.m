clear all; close all; clc;

n = 400; % Number of nodes
k = 7; % Number of clusters
PrSELF = 0.3; % Probability of connection within block
PrF = 0.05; % Scaling correlated probability
PrConst = 0.01; % Probability of connection outside block
VIS = 0; % Turn visualization on/off
[G, Gp, z] = SBM(n,k,PrSELF,PrF,PrConst);

% Visualization
if (VIS)
VIS_SIG = (1.1^k)*(1/k);
e = 2*pi/k;
xy_seed = linspace(0+e,2*pi,k)';
xy = [];
cc = hsv(k);
h1 = figure;
hold all
for i = 1:k
    n_k = sum(z==i);
    s = xy_seed(i);
    [index, ~] = find(z==i);
    for j = 1:n_k
        i_j = index(j);
        xy(i_j,1) = cos(s)+normrnd(0,VIS_SIG);
        xy(i_j,2) = sin(s)+normrnd(0,VIS_SIG);
        plot(xy(i_j,1),xy(i_j,2),'.','color',cc(i,:),'MarkerSize',25)
    end
end
gplot(G,xy,'black-o')
print(h1,'-djpeg','hw5_graph')
end

SIMS = 500;
L = zeros(k,SIMS);
for b = 1:SIMS
    [G, Gp, z] = SBM(n,k,PrSELF,PrF,PrConst);
    cond = zeros(n,1);
    for K = 1:(k)
        cond = cond + z==K;
        Gk = G(cond,cond);
        Gpk = Gp(cond,cond);
        logl = sum(sum(Gk.*log(Gpk) + (1-Gk).*log(1-Gpk)));
        %fprintf('k: %i\tLog Likelihood: %1.4f\n',K,logl)
        L(K,b) = logl;
    end
end

h2 = figure;
plot(-mean(L,2),'black-')
xlabel('number of clusters')
ylabel('log-likelihood')
print(h2,'-djpeg','hw5_logl')

% sample
% n = 100;
% p = 0.21;
% 
% A = rand(n) < p; % directed graph, one block
% 
% likA = p.^A*(1-p).^(1-A);
% loglikA = A*log(p) + (1-A)*log(1-p);
% 
% sumloglikA = sum(loglikA(:));
