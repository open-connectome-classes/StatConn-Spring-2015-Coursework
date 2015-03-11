% sample
n=100;
p=0.21;

A=rand(n)<p;


% log-likelihood for SBM with K=1
likA=p.^A*(1-p).^(1-A);
loglikA=A*log(p)+(1-A)*log(1-p);


prodlikA=prod(likA(:));
sumloglikA=sum(loglikA(:));

exp(sumloglikA)
