clc; clear all; close all;
%Sample, SBM with 4 blocks, undirected, also assumes each cluster has the
%same number of nodes

%----------
%inputs
clusters=3; %True numer of clusters
test_c=7;  %up to how many clusters should we test?
%------------

blocks=clusters^2;
p=rand(1,blocks)*.9+.1; %the p's range from .1 to 1 cuz 0's ln don't get along
n=clusters*30; %total nodes

for c=1:test_c %Find the loss with multiple cluster estimates. Try estimates from 1 to test_c
    for t=1:10   %How many times to sample so we can take an average loss
        
        %sample a graph with the true p
        temp=rand(n/clusters)<p(1);
        for i=2:blocks
            temp=[temp rand(n/clusters)<p(i)];
        end

        A=temp(:,1:n);
        for i=2:clusters
           A=[A;temp(:,((i-1)*n+1:(i*n)))];
        end
        
        %%
        %liklihood function using a estimated number of clusters.
        
        p_est=zeros(1,c^2);
        cluster_est=kmeans(A,c);
        n_block=zeros(1,c^2);
        for row=1:n
            for col=1:n
                if(A(row,col)==1)
                    p_est(cluster_est(row)+(cluster_est(col)-1)*c)=p_est(cluster_est(row)+(cluster_est(col)-1)*c)+1;
                end
                n_block(cluster_est(row)+(cluster_est(col)-1)*c)=n_block(cluster_est(row)+(cluster_est(col)-1)*c)+1;
            end
        end

        p_est=p_est./n_block;
        
        AP=zeros(n,n);
        for row=1:n
            for col=1:n
                    AP(row,col)=p_est(cluster_est(row)+(cluster_est(col)-1)*c);
            end
        end
        
        
        loglikA = A.*log(AP)+(1-A).*log(1-AP);
        sumloglikA=sum(sum(loglikA(:)));
        %loss_4block(t)=exp(sumloglikA);
        loss_cblock(c,t)=sumloglikA;
    end
    %graph the graph estimated with test_c clusters of equal sizes
    figure(c);
    imagesc(AP); colormap(gray);
end

%%
%plot stuff and find risk
figure(c+1)
imagesc(A); colormap(gray);
title('true Graph');

risk=-mean(loss_cblock,2); %risk of all the cluster estimates

figure(c+2);
plot(risk)