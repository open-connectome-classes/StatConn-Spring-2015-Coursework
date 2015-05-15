function Modularity_Cluster_Seizure(preictal_adj_mat,ictal_adj_mat,p)
    % % Uncomment these lines if you would like to plot original
    % % adjacency matrices for each patient
    % figure(1)
    % subplot(1,2,1)
    % matplot(preictal_adj_mat)
    % subplot(1,2,1)
    % matplot(ictal_adj_mat)

    % % Cluster preictal network
    [r,~] = size(ictal_adj_mat);
    COMTY1 = cluster_jl(preictal_adj_mat);
    commID1 = COMTY1.COM{1};
    comms1 = {};
    comms1{size(COMTY1.SIZE{1},2)} = [];

    for i = 1:size(commID1,2)
        comms1{commID1(i)} = [comms1{commID1(i)} i];
    end

    maximum1 = 0;
    for i = 1:size(comms1,2)
        if size(comms1{i},2) > maximum1
            maximum1 = size(comms1{i},2);
        end
    end

    clusters1 = zeros(size(comms1,2),maximum1);
    for i = 1:size(comms1,2)
        clusters1(i,1:size(comms1{i},2)) = comms1{i};
    end

    mod1 = COMTY1.MOD(1);

    % % Cluster ictal network
    COMTY2 = cluster_jl(ictal_adj_mat);
    commID2 = COMTY2.COM{1};
    comms2 = {};
    comms2{size(COMTY2.SIZE{1},2)} = [];

    for i = 1:size(commID2,2)
        comms2{commID2(i)} = [comms2{commID2(i)} i];
    end

    maximum2 = 0;
    for i = 1:size(comms2,2)
        if size(comms2{i},2) > maximum2
            maximum2 = size(comms2{i},2);
        end
    end

    clusters2 = zeros(size(comms2,2),maximum2);
    for i = 1:size(comms2,2)
        clusters2(i,1:size(comms2{i},2)) = comms2{i};
    end

    mod2 = COMTY2.MOD(1);
    
    % % Visualize the clusters
    % Generate reshuffled preictal adjacency matrix
    labels1 = [];
    for i = 1:size(clusters1,1)
        x = clusters1(i,:);
        x(x==0) = [];
        labels1 = [labels1 x];
    end
    
    preictal_shuffled = zeros(size(preictal_adj_mat));
    for i = 1:size(labels1,2)
        r = labels1(i);
        for j = 1:size(labels1,2)
            c = labels1(j);
            preictal_shuffled(i,j) = preictal_adj_mat(r,c);
        end
    end
    
    % Generate reshuffled ictal adjacency matrix
    labels2 = [];
    for i = 1:size(clusters2,1)
        x = clusters2(i,:);
        x(x==0) = [];
        labels2 = [labels2 x];
    end
    
    ictal_shuffled = zeros(size(ictal_adj_mat));
    for i = 1:size(labels2,2)
        r = labels2(i);
        for j = 1:size(labels2,2)
            c = labels2(j);
            ictal_shuffled(i,j) = ictal_adj_mat(r,c);
        end
    end
    
    if p == 1
        subplot(1,2,1)
        matplot(preictal_shuffled)
        title('Preictal')
        subplot(1,2,2)
        title('Ictal')
        matplot(ictal_shuffled)
    end
end