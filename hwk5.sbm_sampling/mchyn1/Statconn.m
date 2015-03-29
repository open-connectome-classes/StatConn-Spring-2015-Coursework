clear all;close all;clc;

%SBM
 %{
%sample SBM of 100 vertices and 2 blocks
n = 100;
p1=0.31;
p2 = 0.69;
p12 = .10;

A = rand(n)<p1;
B = rand(100,50)<p12;
C = rand(50,100)<p12;
D = rand(50,50)<p2;

M = [A B ; C D];

spy(M)
%}
%sample
n = 500; %num verticies
for realc = 2:5 %num clusters
[M Preal] = makeSBM(n,realc);

phat=sum(M(:))/n^2;

%likelihood with 1 cluster
%{
likA=nan(n);
for u=1:n
    for v=1:n
        likA(u,v)= p^A(u,v)*(1-p)^(1-A(u,v));
    end
end
%}
likM=phat.^M*(1-phat).^(1-M);
loglikM=M*log(phat)+(1-M)*log(1-phat);
prodlikM=prod(likM(:));
sumloglikM = nan(3,1);
sumloglikM(1)=sum(loglikM(:));
a = exp(sumloglikM);

%repeat over k blocks
%{
clust = kmeans(M,2);
aa = M(clust==1, clust==1);
bb = M(clust==2, clust==2);
ab = M(clust==1, clust==2);
ba = M(clust==2, clust==1);
M2 = [aa ab;ba bb];

likM2 = nan(150);
loglikM2 = nan(150);
for i = 1:150
    for j = 1:150
        if i <= 100 && j <= 100
            likM2(i,j)=p1.^M2(i,j)*(1-p1).^(1-M2(i,j));
            loglikM2(i,j)=M2(i,j)*log(p1)+(1-M2(i,j))*log(1-p1);
        elseif i>100 && j>100
            likM2(i,j)=p2.^M2(i,j)*(1-p2).^(1-M2(i,j));
            loglikM2(i,j)=M2(i,j)*log(p2)+(1-M2(i,j))*log(1-p2);
        else
            likM2(i,j)=p12.^M2(i,j)*(1-p12).^(1-M2(i,j));
            loglikM2(i,j)=M2(i,j)*log(p12)+(1-M2(i,j))*log(1-p12);
        end
    end
end
prodlikM2=prod(likM2(:));
sumloglikM(2,1)=sum(loglikM2(:));
%}


for c = 1:10
clust = kmeans(M,c);
m= cell(c,c);
for i = 1:c
    for j = 1:c
        m{i,j} = M(clust==i,clust ==j);
    end
end
%{
m{1,1} = M(clust==1, clust==1);
a22 = M(clust==2, clust==2);
a33 = M(clust==3, clust==3);
a12 = M(clust==1, clust==2);
a13 = M(clust==1, clust==3);
a21 = M(clust==2, clust==1);
a23 = M(clust==2, clust==3);
a31 = M(clust==3, clust==1);
a32 = M(clust==3, clust==2);
%}
X = nan(c*c,2);
in = 1;
for j = 1:c
    for k = 1:c
        X(in,:) = size(m{j,k});
        in = in+1;
    end
end
in = 1;
M3 = cell2mat(m);
%figure
%spy(M3)
P = nan(c,c);
for i = 1:c
    for j = 1:c
        P(i,j) = sum(m{i,j}(:))/(X(in,1)*X(in,2));
        in = in+1;
    end
end
loglik = cell(c,c);
X = cell(c,c);
Y = cell(c,c);
for i = 1:c
    for j = 1:c
        X{i,j} = P(i,j).^m{i,j};
        Y{i,j} = (1-P(i,j)).^(1-m{i,j});
        loglik{i,j}=m{i,j}*log(P(i,j))+(1-m{i,j})*log(1-P(i,j));
    end
end
likM3 = cell2mat(X)*cell2mat(Y);
loglikM3 = cell2mat(loglik);
prodlikM3=prod(likM3(:));
sumloglikM(c)=sum(loglikM3(:));
end

figure
plot(1:length(sumloglikM),sumloglikM);
xlabel('cluster size');
ylabel('log likelihood');
title(sprintf('actual num clusters = %d',realc));


end