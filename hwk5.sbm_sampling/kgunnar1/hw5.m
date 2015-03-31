% Connectomics

n = 100;
p = 0.21;

A = rand(n)<p;

%sum(A(:))/n^2;

% Likelihood
% In general we wan't to compute log likelihoods because likelihoods are
% small numbers and when multiply very many very small numbers we get tiny
% numbers and get numerical rounding errors or even run out of...
% So in practice - compute Log Likelihood

% % Log likelihood:
% This is just the likelihood
% likA = nan(n);
% for u = 1:n
%     for v = 1:n
%         likA(u,v) = p^A(u,v)*(1-p)^(1-A(u,v));
%     end
% end

% or don't need loop:
% The likelihood
likA = p.^A*(1-p).^(1-A);

% The log likelihood
loglikA = A*log(p)+(1-A)*log(1-p);

% Product of everything:
prodLikA = prod(likA(:));
sumLogLikA = sum(loglikA(:));
% then probability
probA = exp(sumLogLikA)

% Model 2
n = 100;
p1 = 0.25;
p2 = 0.42;
p3 = 0.12;
p4 = 0.6;

B = zeros(n);
B(1:25,1:25) = rand(25)<p1;
B(26:50,26:50) = rand(25)<p2;
B(51:75,51:75) = rand(25)<p3;
B(76:100,76:100) = rand(25)<p4;
imagesc(B)

% Likelihood and log likelihood 
likB = p1.^B(1:25,1:25)*(1-p1).^(1-B(1:25,1:25))*p2.^B(26:50,26:50)*(1-p2).^(1-B(26:50,26:50))*p3.^B(51:75,51:75)*(1-p3).^(1-B(51:75,51:75))*p4.^B(76:100,76:100)*(1-p4).^(1-B(76:100,76:100))
loglikB = B(1:25,1:25)*log(p1)+ log(1-p1)*(1-B(1:25,1:25))+ log(p2)*B(26:50,26:50)+log(1-p2)*(1-B(26:50,26:50))+log(p3)*B(51:75,51:75)+log(1-p3)*(1-B(51:75,51:75))+log(p4)*B(76:100,76:100)+log(1-p4)*(1-B(76:100,76:100))

%product of likelihood/sum of log likelihood
prodlikB=prod(likB(:));
sumloglikB=sum(loglikB(:));

probB = exp(sumloglikB)


% Assuming only 2 clusters instead of four
temp1 = B(1:50,1:50);
temp2 = B(51:100,51:100);
p1_1 = sum(temp1(:))/50^2;
p2_1 = sum(temp2(:))/50^2;

loglikB_2 = B(1:50,1:50)*log(p1_1)+ log(1-p1_1)*(1-B(1:50,1:50))+ log(p2_1)*B(51:100,51:100)+log(1-p2_1)*(1-B(51:100,51:100));
sumloglikB_2 = sum(loglikB_2(:));

probB_2 = exp(sumloglikB_2)         % <--- Get a smaller likelihood when using the wrong number of clusters