clear all; close all; clc;

fprintf('#----------------------------------------------------------------#\n')
fprintf('| CWAS Statistical Sandbox :: Version 1.0                        |\n')
fprintf('| Jerry Wang (jwang158@jhu.edu)                                  |\n')
fprintf('#----------------------------------------------------------------#\n')
fprintf('Please read the README for more details.\n\n')

n = 1500; % Number of nodes
k = 50; % Number of clusters
NP = 100; % Number of patients
TOP = 20; % Number of top-correlated features to cherry pick
N_PERM = 10000; % Number of permutation iterations
PrSELF = 0.3; % Probability of connection within block
PrF = 0.8; % Scaling correlated probability
PrConst = 0.0001; % Probability of connection outside block
VIS = 0; % Turn visualization on/off
VIS_PERM = 0; % Turn permutation histogram visualization on/off

% Create random phenotypes independent of features
Y = zeros(NP,1);
Y(1:(NP/2),1) = ones(NP/2,1);
Y = Y(randperm(NP));

% Create random features
X = zeros(NP,n);
fprintf('[*] Simulating random data for %i patients at %i voxels each\n',NP,n)
for i = 1:NP
    [~, Gp, ~] = SBM(n,k,PrSELF,PrF,PrConst);
    X(i,:) = sum(Gp);
    
    % Print progress
    if (i == 1)
        fprintf('    Percent: 0..')
    elseif (mod(i-1,NP/10) == 0)
        fprintf('%i..',100*(i-1)/NP)
    end
end
fprintf('\n')

% Run cherry pick correlations
fprintf('[*] Running correlations and picking top %i features...\n',TOP)
R = corr(X,Y);
[~,Ri] = sort(R.^2,'descend');
top_features = Ri(1:TOP);
X2 = X(:,top_features);

% Permutation testing
fprintf('[*] Permutation testing with %i iterations...\n', N_PERM)
Pperm = zeros(TOP,1);
Xperm = X(:,randperm(n));
Xperm = Xperm(:,1:TOP); % Randomly cherry pick TOP number of features
for i = 1:TOP
    % The absolute difference in mean statistic
    m_test = abs(mean(Xperm(Y==1,i)) - mean(Xperm(Y==0,i)));
    parfor j = 1:N_PERM
        % Permute the phenotype labels
        Yperm = Y(randperm(NP));
        M(j) = abs(mean(Xperm(Yperm==1,i)) - mean(Xperm(Yperm==0,i)));
    end
    
    % Calculate the one-sided p-value by counting
    Pperm(i) = sum(M > m_test)/N_PERM;

    % Print progress
    if (i == 1)
        fprintf('    Percent: 0..')
    elseif (mod(i-1,TOP/10) == 0)
        fprintf('%i..',100*(i-1)/TOP)
    end
    
    % Optional visualization of the permutation test
    if (VIS_PERM)
        [a,b] = hist(M);
        figure;
        plot(b,a,'black-o')
        hold on
        plot([m_test m_test],[0 max(a)],'red--')
    end
end
fprintf('\n')

% Perform ANOVA test
fprintf('[*] Performing ANOVA test...\n')
X = X2;
P = zeros(TOP,1);
for i = 1:TOP
    P(i) = anova1([X(Y==0,i),X(Y==1,i)],[0 1],'off');
end

% Print output
fprintf('[!] Done.\n\n')
fprintf('Number of Significant nodes by CWAS: %i\n\n',sum(P<0.05/TOP))
fprintf('Number of Significant nodes by proposed method: %i\n\n',sum(Pperm<0.05/TOP))

% Print graphs
h_1 = figure;
subplot(2,1,1) % Upper subplot
plot(-log10(P),'blacko')
hold on
plot([1 TOP],[-log10(0.05/TOP),-log10(0.05/TOP)],'red--')
title('Significant features found by CWAS (FWER < 0.05)')
ylabel('-log10(p-value)')
xlabel('Top-correlated feature')
legend('p-value from ANOVA','Bonferroni significance threshold')
subplot(2,1,2) % Lower subplot
plot(-log10(sort(Pperm)+1e-6),'blacko')
hold on
plot([1 TOP],[-log10(0.05/TOP),-log10(0.05/TOP)],'red--')
title('Significant features found by proposed method (FWER < 0.05)')
ylabel('-log10(p-value)')
xlabel('Top-correlated feature')
legend('p-value from permutation test','Bonferroni significance threshold')
% Write to disk
print(h_1,'-djpeg','p-values_CWAS')

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
end

fprintf('\n--- EXIT SUCCESS ---\n')