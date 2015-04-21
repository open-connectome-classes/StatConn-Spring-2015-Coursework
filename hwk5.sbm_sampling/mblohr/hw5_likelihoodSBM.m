% calculate log-likelihood of SBM

% sample
n = 16; 
n11 = 8;
n12 = 8;
n21 = 8;
n22 = 8
p = [0.4 0.2; 0.2 0.4];
pMatrix = [p(1,1)*ones(n11) p(1,2)*ones(n12);  ...
            p(2,1)*ones(n21) p(2,2)*ones(n22)]
A = rand(n)<pMatrix

% log-likelihood for SBM
likA = (pMatrix.^A).*(ones(n)-pMatrix).^(ones(n)-A)
loglikA = A.*log(pMatrix) + (ones(n)-A).*log(ones(n)-pMatrix)

