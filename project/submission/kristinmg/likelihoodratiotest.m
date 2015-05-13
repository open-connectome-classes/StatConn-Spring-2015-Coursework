function [confusion_matrix, accuracy] = likelihoodratiotest(data, num1, num2);

%data - for example degree of a node
% num1 - number of patients in group1
% num2 - number of patients in group2

num_bins = 10;      

patient_group = zeros(1,length(data));  % The classification of each individual
p1_all = nan(1,length(data));           % Probability of being a Female
p2_all = nan(1,length(data));           % Probability of being a Male

for patient = 1:length(data)
    % Leave the individual we are classifying out
    patients_used = data(1:end ~=patient);
    
    if patient <=num1     % If patient left out is a Female
        % Calculate priors for each group
        prior1 = (num1-1)/(num1+num2-1);        % Prior for Females.
        prior2 = num2/(num1+num2-1);            % Prior for Males
        
        % Compute distribution for each group
        group1 = patients_used(1:num1-1);           % Left one Female out so have (num1-1)Females
        group2 = patients_used(num1:end);         % num2 Males
        [N1, edges1] = histcounts(group1,num_bins,'Normalization','pdf');       % Compute distribution of Females
        [N2, edges2] = histcounts(group2,num_bins,'Normalization','pdf');       % Compute distribution of Males
        
        % Find the center of each bar in the distributions
        for c1 = 2:length(edges1)
            centers1(c1-1) = edges1(c1-1)+(edges1(c1)-edges1(c1-1))/2;
        end
        
        for c2 = 2:length(edges2)
            centers2(c2-1) = edges2(c2-1)+(edges2(c2)-edges2(c2-1))/2;
        end
        
        patient_observed = data(patient);   % The value for the patient I'm classifying
        
        % Probability of patient being in group1 (Females)
            [~, index1] = min(abs(patient_observed-centers1));      % Find the index of the bin closest to the value of our patient being classified
            if index1 == 1 % If the closest bin is one of the end bins we check if the value of our patient falls within the distribution
                if patient_observed < edges1(1)     % If patient falls outside the distribution (smaller value than the leftmost edge)
                    p1 = 0;
                else
                    p1 = N1(index1);
                end
            elseif index1 == length(N1)     % If the closest bin is one of the end bins we check if the value of our patient falls within the distribution             
                if patient_observed > edges1(end)         % If patient falls outside the distribution (higher value than the rightmost edge)
                p1 = 0;
                 else
                    p1 = N1(index1);
                end
            else      
                % For some bins the probability is 0 in the histogram,
                % then finding the mean of the two closest probabilities on
                % each side of the zero bin
                if N1(index1) == 0 && index1 > 1
                    p1 = abs(N1(index1+1)+N1(index1-1))/2;
                else
                    p1 = N1(index1);
                end
            end
        
            % Probability of patient being in group2 (Males)
            [~, index2] = min(abs(patient_observed-centers2));          
            if index2 == 1  % If the closest bin is one of the end bins we check if the value of our patient falls within the distribution
                if patient_observed < edges2(1) % If patient falls outside the distribution (smaller value than the leftmost edge)
                    p2 = 0;
                else
                    p2 = N2(index2);
                end
            elseif index2 == length(N2)        % If the closest bin is one of the end bins we check if the value of our patient falls within the distribution          
                if patient_observed > edges2(end)      % If patient falls outside the distribution (higher value than the rightmost edge)    %
                p2 = 0;
                 else
                    p2 = N2(index2);
                end
            else 
                % For some bins the probability is 0 in the histogram,
                % then finding the mean of the two closest probabilities
                if N2(index2) == 0 && index2 > 1
                    p2 = abs(N2(index2+1)+N2(index2-1))/2;
                else
                    p2 = N2(index2);
                end
            end
            
    else        % If patient being classified is in group2 (Males)
        % Calculate priors for each group
        prior1 = num1/(num1+num2-1);        % Now we have num1 Females
        prior2 = (num2-1)/(num1+num2-1);    % And num2-1 MAles
        
        % Compute distribution for each group
        group1 = patients_used(1:num1);
        group2 = patients_used(num1+1:end);
        [N1, edges1] = histcounts(group1,num_bins,'Normalization','pdf');
        [N2, edges2] = histcounts(group2,num_bins,'Normalization','pdf');
        
        for c1 = 2:length(edges1)
            centers1(c1-1) = edges1(c1-1)+(edges1(c1)-edges1(c1-1))/2;
        end
        
        for c2 = 2:length(edges2)
            centers2(c2-1) = edges2(c2-1)+(edges2(c2)-edges2(c2-1))/2;
        end
        
        patient_observed = data(patient);
        % Probability of patient being in group1 (Females)
            [~, index1] = min(abs(patient_observed-centers1));
           if index1 == 1  % If the closest bin is one of the end bins we check if the value of our patient falls within the distribution
                if patient_observed < edges1(1) % If patient falls outside the distribution (lower value than the leftmost edge)    %
                    p1 = 0;
                else
                    p1 = N1(index1);
                end
            elseif index1 == length(N1)                  % If the closest bin is one of the end bins we check if the value of our patient falls within the distribution
                if patient_observed > edges1(end)         % If patient falls outside the distribution (higher value than the rightmost edge)    %
                p1 = 0;
                 else
                    p1 = N1(index1);
                end
            else      
                % For some bins the probability is 0 in the histogram,
                % then finding the mean of the two closest probabilities
                if N1(index1) == 0 && index1 > 1
                    p1 = abs(N1(index1+1)+N1(index1-1))/2;
                else
                    p1 = N1(index1);
                end
            end
            
            [~, index2] = min(abs(patient_observed-centers2));           
            if index2 == 1    % If the closest bin is one of the end bins we check if the value of our patient falls within the distribution
                if patient_observed < edges2(1)             % If patient falls outside the distribution (lower value than the leftmost edge)    %
                    p2 = 0;
                else
                    p2 = N2(index2);
                end
            elseif index2 == length(N2)                  % If the closest bin is one of the end bins we check if the value of our patient falls within the distribution
                if patient_observed > edges2(end)         % If patient falls outside the distribution (higher value than the rightmost edge)    %
                p2 = 0;
                 else
                    p2 = N2(index2);
                end
            else 
                % For some bins the probability is 0 in the histogram,
                % then finding the mean of the two closest probabilities
                if N2(index2) == 0 && index2 > 1
                    p2 = abs(N2(index2+1)+N2(index2-1))/2;
                else
                    p2 = N2(index2);
                end
            end
    end
    
    % Computing the posterior probabilities of the patient being a Female or Male (likelihood*prior)
    p1_all(patient) = p1*prior1;            % Probability of being in Females group
    p2_all(patient) = p2*prior2;            % Probabiltiy of being in Males group
    
    if p1_all(patient) > p2_all(patient)
        patient_group(patient) = 1;     % Female
    elseif p2_all(patient) > p1_all(patient)
        patient_group(patient) = 2;     % Male
    else % if p1 = p2 = 0 I'll "flip a coin" to classify the patient
        flipcoin = rand();
        if flipcoin <= 0.5
            patient_group(patient) = 1;     % Female
        else
            patient_group(patient) = 2;     % Male
        end
    end
end

% Constructing the confusion matrix
malemale = 0;
malefemale = 0;
femalemale  = 0;
femalefemale = 0;

for i = 1:length(patient_group)    
    if i <= num1  % For true Females
        if patient_group(i) == 1            % Counting all correct Females
            femalefemale = femalefemale+1;
        else
            femalemale = femalemale+1;              % Counting all Females classified as Males
        end
    else     % For true Males
        if patient_group(i) == 2    
            malemale = malemale+1;              % Counting all correct Males
        else
            malefemale = malefemale+1;              % Counting all Males classified as Females
        end
    end
end

% Filling in confusing matrix
confusion_matrix = zeros(2);
confusion_matrix(1,1) = malemale;
confusion_matrix(1,2) = malefemale;
confusion_matrix(2,1) = femalemale;
confusion_matrix(2,2) = femalefemale;

% Calculating the accuracy of the classification
accuracy = trace(confusion_matrix)/length(patient_group);

% % Spitting out the patients that have p1 = p2 = 0
% for i = 1:length(p1_all)
% if p1_all(i) == 0 && p2_all(i) == 0
% fprintf('For patient #%d, p1 = p2 = %.2f\n',i,p1_all(i)) 
% end
% end
        