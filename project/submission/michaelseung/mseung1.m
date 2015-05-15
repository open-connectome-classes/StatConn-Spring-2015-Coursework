% Michael Seung %
% Statistical Connectomics %
% Final Project %

clear all;
url='https://www.dropbox.com/sh/idt3d0gylplyo31/AAAecrQu5UV0DV8EXbnbI64Ia/mseung.zip?dl=0';
cmd=['wget ' url ' -O "mseung.zip" --no-check-certificate'];
system(cmd);
unzip('mseung.zip');
addpath('mseung');


fig=importdata('dat.mat');
fig_mixed=importdata('dat2.mat');

figure(1)
imagesc(fig)
figure(2)
imagesc(fig_mixed)

[In1,O1] = fmincon(@(p) -sum(sum(fig*log(p)+(1-fig)*log(1-p))),0.01,[],[],[],[],0,1);

% SBM with 2 clusters
aa = kmeans(fig,2);
ab = fig(aa==1,aa==1);
ac = fig(aa==2,aa==2);
ad = fig(aa==1,aa==2);
[In2,O2] = fmincon(@(p) -sum(sum(ab*log(p(1))+(1-ab)*log(1-p(1))))-...
    sum(sum(ac*log(p(2))+(1-ac)*log(1-p(2))))-...
    sum(sum(ad*log(p(3))+(1-ad)*log(1-p(3))))...
    ,[0.01 0.01 0.001],[],[],[],[],[0 0 0],[1 1 1])

% /w 3 clusters
ba = kmeans(fig,3);
bb = fig(ba==1,ba==1);
bc = fig(ba==2,ba==2);
bd = fig(ba==3,ba==3);
be = fig(ba==1,ba==2);
bf = fig(ba==1,ba==3);
bg = fig(ba==2,ba==3);
[In3,O3] = fmincon(@(p) -sum(sum(bb*log(p(1))+(1-bb)*log(1-p(1))))-...
    sum(sum(bc*log(p(2))+(1-bc)*log(1-p(2))))-...
    sum(sum(bd*log(p(3))+(1-bd)*log(1-p(3))))-...
    sum(sum(be*log(p(4))+(1-be)*log(1-p(4))))-...
    sum(sum(bf*log(p(5))+(1-bf)*log(1-p(5))))-...
    sum(sum(bg*log(p(6))+(1-bg)*log(1-p(6))))...
    ,[0.01 0.01 0.01 0.001 0.001 0.001],[],[],[],[],[0 0 0 0 0 0],[1 1 1 1 1 1])


[In1_mixed,O1_mixed] = fmincon(@(p) -sum(sum(fig_mixed*log(p)+(1-fig_mixed)*log(1-p))),0.01,[],[],[],[],0,1);

% SBM with 2 clusters
aa = kmeans(fig_mixed,2);
ab = fig_mixed(aa==1,aa==1);
ac = fig_mixed(aa==2,aa==2);
ad = fig_mixed(aa==1,aa==2);
[In2_mixed,O2_mixed] = fmincon(@(p) -sum(sum(ab*log(p(1))+(1-ab)*log(1-p(1))))-...
    sum(sum(ac*log(p(2))+(1-ac)*log(1-p(2))))-...
    sum(sum(ad*log(p(3))+(1-ad)*log(1-p(3))))...
    ,[0.01 0.01 0.001],[],[],[],[],[0 0 0],[1 1 1])

% /w 3 clusters
ba = kmeans(fig_mixed,3);
bb = fig_mixed(ba==1,ba==1);
bc = fig_mixed(ba==2,ba==2);
bd = fig_mixed(ba==3,ba==3);
be = fig_mixed(ba==1,ba==2);
bf = fig_mixed(ba==1,ba==3);
bg = fig_mixed(ba==2,ba==3);
[In3_mixed,O3_mixed] = fmincon(@(p) -sum(sum(bb*log(p(1))+(1-bb)*log(1-p(1))))-...
    sum(sum(bc*log(p(2))+(1-bc)*log(1-p(2))))-...
    sum(sum(bd*log(p(3))+(1-bd)*log(1-p(3))))-...
    sum(sum(be*log(p(4))+(1-be)*log(1-p(4))))-...
    sum(sum(bf*log(p(5))+(1-bf)*log(1-p(5))))-...
    sum(sum(bg*log(p(6))+(1-bg)*log(1-p(6))))...
    ,[0.01 0.01 0.01 0.001 0.001 0.001],[],[],[],[],[0 0 0 0 0 0],[1 1 1 1 1 1])

% Figure
figure(3)
plot(1:3,-[O1_mixed O2_mixed O3_mixed],'.r', 'MarkerSize',40);
hold all;
plot(1:3,-[O1 O2 O3],'.b', 'MarkerSize',40);
xlabel('Number of Clusters');
ylabel('Log Likelihood');
legend('Mixed Model','Not Mixed','Location','northwest')


