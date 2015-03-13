n_a = 500;
n_b = 500;
p_aa = 0.05;
p_ab = 0.01;
p_bb = 0.05;

y = zeros(n_a+n_b);

% SBM
y(1:n_a,1:n_a) = rand(n_a)<p_aa;
y(n_a+1:n_b+n_a,n_a+1:n_b+n_a) = rand(n_b)<p_bb;
y(1:n_a,n_a+1:n_b+n_a) = rand(n_a,n_b)<p_ab;
y(n_a+1:n_b+n_a,1:n_a) = y(1:n_a,n_a+1:n_b+n_a);


% 1 cluster
[lhood_1,val_1] = fmincon(@(p) -sum(sum(y*log(p)+(1-y)*log(1-p))),0.01,...
    [],[],[],[],0,1)

% 2 cluster
y_2 = kmeans(y,2);
b_11 = y(y_2==1,y_2==1);
b_22 = y(y_2==2,y_2==2);
b_12 = y(y_2==1,y_2==2);
[lhood_2,val_2] = fmincon(@(p) -sum(sum(b_11*log(p(1))+(1-b_11)*log(1-p(1))))-...
    sum(sum(b_22*log(p(2))+(1-b_22)*log(1-p(2))))-...
    sum(sum(b_12*log(p(3))+(1-b_12)*log(1-p(3))))...
    ,[0.01 0.01 0.001],[],[],[],[],[0 0 0],[1 1 1])

% 3 clusters
y_3 = kmeans(y,3);
c_11 = y(y_3==1,y_3==1);
c_22 = y(y_3==2,y_3==2);
c_33 = y(y_3==3,y_3==3);
c_12 = y(y_3==1,y_3==2);
c_13 = y(y_3==1,y_3==3);
c_23 = y(y_3==2,y_3==3);
[lhood_3,val_3] = fmincon(@(p) -sum(sum(c_11*log(p(1))+(1-c_11)*log(1-p(1))))-...
    sum(sum(c_22*log(p(2))+(1-c_22)*log(1-p(2))))-...
    sum(sum(c_33*log(p(3))+(1-c_33)*log(1-p(3))))-...
    sum(sum(c_12*log(p(4))+(1-c_12)*log(1-p(4))))-...
    sum(sum(c_13*log(p(5))+(1-c_13)*log(1-p(5))))-...
    sum(sum(c_23*log(p(6))+(1-c_23)*log(1-p(6))))...
    ,[0.01 0.01 0.01 0.001 0.001 0.001],[],[],[],[],[0 0 0 0 0 0],[1 1 1 1 1 1])

% Figure 1
figure(1)
imagesc(y);
colormap(gray);
title('Adjacency Matrix');

% Figure 2
figure(2)
plot(1:3,-[val_1 val_2 val_3], 'o');
xlabel('Number of Clusters');
ylabel('Log of Likelihood');
colormap(gray);