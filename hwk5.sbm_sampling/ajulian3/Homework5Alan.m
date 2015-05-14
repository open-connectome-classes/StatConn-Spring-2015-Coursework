% Connectomics- Homework 5



% First we are creating a sample. Considering neuronal orientation is
% split into 4 quadrants, let p=.25
x = 100;
p = 0.25;

G = rand(x)<p;
% We must calculate the likelihood. The likelihood function that is best
% for this problem is loglikelihood. Loglikelihood makes the process
% computationally convenient and makes the data more avaialable to
% interpretation. Loglikelihood is practical when examining lower p values.

likA = p.^G*(1-p).^(1-G);

% The log likelihood function 
loglikA = G*log(p)+(1-G)*log(1-p);

% Product of everything:
prodLikA = prod(likA(:)); %since its not log -> summation
sumLogLikA = sum(loglikA(:)); % log -> product
% To determine probability from the log likelihood function, we must
% exponate the sum
probA = exp(sumLogLikA)

% Model 2
n = 100;

% Presented various probabilities that reflect multiple p values
p1 = 0.27;
p2 = 0.47;
p3 = 0.16;
p4 = 0.65;

B = zeros(n);
B(1:25,1:25) = rand(25)<p1;
B(26:50,26:50) = rand(25)<p2;
B(51:75,51:75) = rand(25)<p3;
B(76:100,76:100) = rand(25)<p4;
imagesc(B)

% Likelihood and log likelihood 
likB = p1.^B(1:25,1:25)*(1-p1).^(1-B(1:25,1:25))*p2.^B(26:50,26:50)*(1-p2).^(1-B(26:50,26:50))*p3.^B(51:75,51:75)*(1-p3).^(1-B(51:75,51:75))*p4.^B(76:100,76:100)*(1-p4).^(1-B(76:100,76:100))
loglikB = B(1:25,1:25)*log(p1)+ log(1-p1)*(1-B(1:25,1:25))+ log(p2)*B(26:50,26:50)+log(1-p2)*(1-B(26:50,26:50))+log(p3)*B(51:75,51:75)+log(1-p3)*(1-B(51:75,51:75))+log(p4)*B(76:100,76:100)+log(1-p4)*(1-B(76:100,76:100))

%product of likelihood

prodlikB=prod(likB(:));

% Sum of likelihood
sumloglikB=sum(loglikB(:));

% Have to exponate to determine probability 
probB = exp(sumloglikB)


% Instead of using the four clusters used previously, now we are going to
% run the code with using 2 clusters instead of the 4

% Cluster 1
C1 = B(1:50,1:50);

%Cluster 2
C2 = B(51:100,51:100);

% Defining new probabilities


% New P1 
p1new = sum(C1(:))/50^2;


% New P2
p2new = sum(C2(:))/50^2;

loglikBnew = B(1:50,1:50)*log(p1new)+ log(1-p1new)*(1-B(1:50,1:50))+ log(p2new)*B(51:100,51:100)+log(1-p2new)*(1-B(51:100,51:100));
sumloglikBnew = sum(loglikBnew(:));


probBnew = exp(sumloglikBnew)        
% Interestingly, when the incorrect number of clusters is implemented, we
% notice that the likelihood decreases. 