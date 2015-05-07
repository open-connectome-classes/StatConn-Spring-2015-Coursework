function [ G, Gp, z ] = SBM( n, k, P_SELF, F, Const )
% SBM - Stochastic Block Model
%   Returns a graph G from the SBM model with size n and k clusters
%   Jiarui Wang :: mr.jiarui.wang@gmail.com

    % Constants
    %P_SELF = 0.5; % Probability of connecting with self
    %F = 0.9; % Fraction scaled from self

    % Create root matrices
    z = randi([1 k],n,1);
    Gp = zeros(n,n);
    M = (P_SELF-Const)*eye(k);
    scale = P_SELF;
    for c = 1:(k-1)
        scale = scale * F;
        M = M + diag(scale*ones(k-c,1),c) + diag(scale*ones(k-c,1),-c);
    end
    M = M + Const*ones(k,k);

    % Pull probabilities out of the M root matrix
    parfor i = 1:n;
        for j = 1:n
            Gp(i,j) = M(z(i),z(j));
        end
    end
    G = rand(n) <= Gp;
end

