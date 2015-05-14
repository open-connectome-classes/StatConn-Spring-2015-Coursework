%Sandya Subramanian
%Statistical Connectomics Final Project 
%Spring 2015

%Main Script
%% Load data

loadData();
load('networks.mat');
load('event_labels.mat');
%% Calculate parameters of SBM for each event and plot the difference + calculate accuracy

%SBM based on RR and not RR 
%For each seizure, move the RR block to the top left, binarize edges by thresholding, and compute the
%p(edge|RR) and p(edge|not RR) based on edge density (MLE estimate)
SBM_params = zeros(2,length(event_labels));
%Threshold for edges
edge_thresh = 0.6;
for i = 1:length(event_labels)
    %Get the RR_electrode labels
    RR_elec = networks{2,i};
    %Get the network
    network = networks{1,i};
    %Get the number of electrodes
    num_elec = size(network,1);
    %Get the non-RR electrodes
    other_elec = setdiff(1:num_elec,RR_elec);
    
    %Move the RR electrodes to the top left in the graph
    shuffled_network = network([RR_elec,other_elec],[RR_elec,other_elec]);
    %Threshold the network
    thresh_network = (shuffled_network > edge_thresh);
    %Compute p(edge|RR)
    SBM_params(1,i) = sum(sum(thresh_network(1:length(RR_elec),1:length(RR_elec))))/(length(RR_elec)^2);
    %Compute p(edge|not RR)
    SBM_params(2,i) = sum(sum(thresh_network(end-length(other_elec)+1:end,end-length(other_elec)+1:end)))/(length(other_elec)^2);
    
end

diff_params = SBM_params(1,:) - SBM_params(2,:);

%Plot difference: p(edge|RR) - p(edge|not RR)
figure;
gscatter(1:length(SBM_params),diff_params,event_labels);
hold on;
legend('Failures','Successes')
xlabel('Seizure Number')
ylabel('P(Edge|RR) - P(Edge|Not RR)')
title('Difference in SBM Probabilities for Each Seizure')

%% Leave one out CV
%Each time compute the threshold that works best

predictions = zeros(size(event_labels));
for n = 1:length(event_labels)
    %Leave out the nth entry in both diff_params and event_labels
    train_set = diff_params(setdiff(1:length(event_labels),n));
    train_labels = event_labels(setdiff(1:length(event_labels),n));
    
    %Compute the best threshold between -1 and 1
    best_t = -1;
    best_acc = 0;
    for t = -1:0.01:1
        check = train_set < t;
        acc = sum(check' == train_labels)/length(train_labels);
        if acc > best_acc
            %If better accuracy, save the threshold
            best_t = t;
            best_acc = acc;
        end
    end
       
    %Predict the last outcome
    predictions(n) = (diff_params(n) < best_t);
end

%% Compute results: accuracy

%Compute accuracy
SBM_accuracy = sum(predictions == event_labels)/length(event_labels)

%% Run random assignment test and check accuracy

%For 10000 iterations - randomly assign 14 successes and 25 failures and
%check the accuracy

num_iter = 10000;
track_accuracy = zeros(1,num_iter);
num_successes = sum(event_labels);

for i = 1:num_iter
    labels = zeros(length(event_labels),1);
    
    %Randomly pull 14 values from 1:39
    sample = randsample(length(event_labels),num_successes);
    
    %Set those to 1
    labels(sample) = 1;
    
    %Compute accuracy
    track_accuracy(i) = sum(labels == event_labels)/length(event_labels);
end

%Get average accuracy across trials
avg_perm_accuracy = mean(track_accuracy)
