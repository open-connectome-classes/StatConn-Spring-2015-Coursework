%% Statistical Connectomics
% Final Project
% Kristin Maria Gunnarsdottir

% Please note - this takes a while to run (around 5 min on my laptop)
%% Get the data 
foldername = 'kristinmg_data.zip';
data_address = ' https://www.dropbox.com/sh/idt3d0gylplyo31/AAD9tcrXsIvE3f5pYPNMRsw2a/kristinmg_data.zip?dl=0';
opts = ' --no-check-certificate';
cmd = strcat('wget', ' -O kristinmg_data.zip', data_address, opts);
system(cmd);
unzip(foldername);
currentdir = pwd;
newdir = strcat(currentdir,'\kristinmg_data');
cd(newdir)

% Get subject data
fileID = fopen('kki42_subjectinformation.csv','r');

headers = textscan(fileID, '%s %s %s %s %s %s %s', 1, 'Delimiter', ',');  % Get headers
data = textscan(fileID, '%f %f %f %s %f %s %s','Delimiter', ',');

fclose(fileID);

subject_num = data{1};
subject_gender = data{6};
female_indices = find(strcmp(subject_gender,'F'));   % 20 Female scans, 22 Male scans
male_indices = find(strcmp(subject_gender,'M'));   % 20 Female scans, 22 Male scans

% Create matrices to store measures
degree = zeros(42,70);
w_degree = degree;
eigen_centrality = degree;
cosine_similarity = zeros(42,70,70);

% For each patient
for patientnum = 1:42
    cd(newdir)
    dataname = strcat('KKI-21_KKI2009-',num2str(patientnum),'_small_graph.mat');
    temp = open(dataname);
    % Get adjacency matrix
    temp_graph = temp.fibergraph;

    % % Visualize adjacency matrix
    % imagesc(temp_graph)

    % Creating a normal adjacency matrix (filling up with zeros)
    A = zeros(size(temp_graph)); 
    [i,j,v] = find(temp_graph);
    for k = 1:length(i)
        A(i(k),j(k))=v(k);
        A(j(k),i(k)) = v(k);        % since graph is undirected, A_ij = A_ji
    end
    
    % Normalize values to lie between 0 and 1
    A = A/max(max(A));
    
    tempdir = strsplit(newdir,'kristinmg_data');
    newdir2 = tempdir{1};
    cd(newdir2)
    degree(patientnum,:) = degree_centrality(A);
    w_degree(patientnum,:) = weighted_degree(A);
    eigen_centrality(patientnum,:) = eigenvector_centrality(A);
    cosine_similarity(patientnum,:,:) = cossim(A);
end

%% Bootstraping to choose features with most significant difference

% Degree
female_degree = degree(female_indices,:);
male_degree = degree(male_indices,:);

pvalues_degree = zeros(size(degree,2),1);
for i = 1:size(degree,2)
    pvalues_degree(i) = my_bootstrap(female_degree(:,i),male_degree(:,i),20,22,10000);
end

% Observe if any p-values < 0.05, i.e. there is a significant difference
% between the two groups

index_significant_degree = find(pvalues_degree < 0.05);           % For degree I want to use nodes 6, 20, 44 as features for likelihood ratio test. Node 6 gives lowest value
p_significant_degree = pvalues_degree(index_significant_degree);

% Weighted Degree
female_w_degree = w_degree(female_indices,:);
male_w_degree = w_degree(male_indices,:);

pvalues_w_degree = zeros(size(w_degree,2),1);
for i = 1:size(w_degree,2)
    pvalues_w_degree(i) = my_bootstrap(female_w_degree(:,i),male_w_degree(:,i),20,22,10000);
end

% Observe if any p-values < 0.05, i.e. there is a significant difference
% between the two groups

index_significant_w_degree = find(pvalues_w_degree < 0.05);           % For weighted degree I want to use nodes 25, 39, 60
p_significant_w_degree = pvalues_w_degree(index_significant_w_degree);

% Eigenvector centrality
female_eig_cent = eigen_centrality(female_indices,:);
male_eig_cent = eigen_centrality(male_indices,:);

pvalues_eig_cent = zeros(size(eigen_centrality,2),1);
for i = 1:size(eigen_centrality,2)
    pvalues_eig_cent(i) = my_bootstrap(female_eig_cent(:,i),male_eig_cent(:,i),20,22,10000);
end

% Observe if any p-values < 0.05, i.e. there is a significant difference
% between the two groups

index_significant_eig_cent = find(pvalues_eig_cent < 0.05);           % For eigenvector centrality I want to use nodes 23, 25, 36
p_significant_eig_cent = pvalues_eig_cent(index_significant_eig_cent);

% Cosine similarity
female_cossim = cosine_similarity(female_indices,:,:);
male_cossim = cosine_similarity(male_indices,:,:);

pvalues_cossim = zeros(size(cosine_similarity,2),size(cosine_similarity,3));
for i = 1:size(cosine_similarity,2)
    pvalues_cossim(i,:) = my_bootstrap(squeeze(female_cossim(:,i,:)),squeeze(male_cossim(:,i,:)),20,22,10000);
end

% Observe if any p-values < 0.05, i.e. there is a significant difference
% between the two groups

[row_significant_cossim, col_significant_cossim] = find(pvalues_cossim  < 0.05);           % For eigenvector centrality I want to use nodes 23, 25, 36
temp1 = row_significant_cossim;
temp2 = col_significant_cossim;

% Since pvalues matrix is symmetric each pair is listed twice. Only want to
% store each pair of nodes once
for i = length(row_significant_cossim):-1:1
    for j = length(col_significant_cossim):-1:1
        if row_significant_cossim(i) == col_significant_cossim(j) && row_significant_cossim(i) ~= row_significant_cossim(j)
            if row_significant_cossim(j) == col_significant_cossim(i)
                row_significant_cossim(j) = [];
                col_significant_cossim(j) = [];
                break
            end
        end
    end
end

indices_significant_cossim = [row_significant_cossim col_significant_cossim];
p_significant_cossim = zeros(1,length(row_significant_cossim));
for i = 1:length(row_significant_cossim)
    p_significant_cossim(i)  = pvalues_cossim(row_significant_cossim(i), col_significant_cossim(i));
end

% Get so many significant for the cosine similarity. Choosing the lowest
% p-values to look at

indices_cossim2 = find(p_significant_cossim < 0.01);

%% Classifying using likelihood ratio test
num1 = 20;  % females
num2 = 22;  % males

% Sorting values so first all females, then males

f_degree = degree(female_indices,:);
m_degree = degree(male_indices,:);
degree = [f_degree; m_degree];

f_w_degree = w_degree(female_indices,:);
m_w_degree = w_degree(male_indices,:);
w_degree = [f_w_degree; m_w_degree];

f_eigcen = eigen_centrality(female_indices,:);
m_eigcen = eigen_centrality(male_indices,:);
eigen_centrality = [f_eigcen; m_eigcen];

f_cossim = cosine_similarity(female_indices,:,:);
m_cossim = cosine_similarity(male_indices,:,:);
cosine_similarity = [f_cossim; m_cossim];

% Constructing confusion matrices and accuracy for all 500 runs
conmat_deg6 = zeros(500,2,2);
accuracy_deg6 = zeros(500,1);
conmat_deg20 = zeros(500,2,2);
accuracy_deg20 = zeros(500,1);
conmat_deg44 = zeros(500,2,2);
accuracy_deg44 = zeros(500,1);
conmat_degall = zeros(500,2,2);
accuracy_degall = zeros(500,1);
conmat_wdeg25 = zeros(500,2,2);
accuracy_wdeg25 = zeros(500,1);
conmat_wdeg39 = zeros(500,2,2);
accuracy_wdeg39 = zeros(500,1);
conmat_wdeg60 = zeros(500,2,2);
accuracy_wdeg60 = zeros(500,1);
conmat_wdegall = zeros(500,2,2);
accuracy_wdegall = zeros(500,1);
conmat_eigcen23 = zeros(500,2,2);
accuracy_eigcen23 = zeros(500,1);
conmat_eigcen25 = zeros(500,2,2);
accuracy_eigcen25 = zeros(500,1);
conmat_eigcen36 = zeros(500,2,2);
accuracy_eigcen36 = zeros(500,1);
conmat_eigcenall = zeros(500,2,2);
accuracy_eigcenall = zeros(500,1);
conmat_cossim6870 = zeros(500,2,2);
accuracy_cossim6870 = zeros(500,1);
conmat_cossim4969 = zeros(500,2,2);
accuracy_cossim4969 = zeros(500,1);
conmat_cossim6668 = zeros(500,2,2);
accuracy_cossim6668 = zeros(500,1);
conmat_cossimall = zeros(500,2,2);
accuracy_cossimall = zeros(500,1);

% Running each test 500 times and finding average accuracy
for i = 1:500
    % Feature: Degree of node 6
    deg6 = degree(:,6);
    [conmat_deg6(i,:,:), accuracy_deg6(i)] = likelihoodratiotest(deg6, num1, num2);

    % Feature: Degree of node 20
    deg20 = degree(:,20);
    [conmat_deg20(i,:,:), accuracy_deg20(i)] = likelihoodratiotest(deg20, num1, num2);

    % Feature: Degree of node 44
    deg44 = degree(:,44);
    [conmat_deg44(i,:,:), accuracy_deg44(i)] = likelihoodratiotest(deg44, num1, num2);

    % Feature: Mean degree of nodes 6,20,44
    deg_all = mean(degree(:,[6,20,44]),2);
    [conmat_degall(i,:,:), accuracy_degall(i)] = likelihoodratiotest(deg_all, num1, num2);
    
    % Feature: Weighted degree of node 25
    w_deg25 = w_degree(:,25);
    [conmat_wdeg25(i,:,:), accuracy_wdeg25(i)] = likelihoodratiotest(w_deg25, num1, num2);

    % Feature: Weighted degree of node 39
    w_deg39 = w_degree(:,39);
    [conmat_wdeg39(i,:,:), accuracy_wdeg39(i)] = likelihoodratiotest(w_deg39, num1, num2);

    % Feature: Weighted degree of node 60
    w_deg60 = w_degree(:,60);
    [conmat_wdeg60(i,:,:), accuracy_wdeg60(i)] = likelihoodratiotest(w_deg60, num1, num2);
    
    % Feature: Mean weighted degree of nodes 25,39,60
    wdeg_all = mean(w_degree(:,[25,39,60]),2);
    [conmat_wdegall(i,:,:), accuracy_wdegall(i)] = likelihoodratiotest(wdeg_all, num1, num2);

    % Feature: Eigenvector centrality of node 23
    eigcen23 = eigen_centrality(:,23);
    [conmat_eigcen23(i,:,:), accuracy_eigcen23(i)] = likelihoodratiotest(eigcen23, num1, num2);

    % Feature: Eigenvector centrality of node 25
    eigcen25 = eigen_centrality(:,25);
    [conmat_eigcen25(i,:,:), accuracy_eigcen25(i)] = likelihoodratiotest(eigcen25, num1, num2);

    % Feature: Eigenvector centrality of node 36
    eigcen36 = eigen_centrality(:,36);
    [conmat_eigcen36(i,:,:), accuracy_eigcen36(i)] = likelihoodratiotest(eigcen36, num1, num2);
    
    % Feature: Mean eigenvector centrality of nodes 23,25,36
    eigcen_all = mean(eigen_centrality(:,[23,25,36]),2);
    [conmat_eigcenall(i,:,:), accuracy_eigcenall(i)] = likelihoodratiotest(eigcen_all, num1, num2);
    

    % Feature: cosine similarity of nodes 68 and 70
    cossim6870 = cosine_similarity(:,68,70);
    [conmat_cossim6870(i,:,:), accuracy_cossim6870(i)] = likelihoodratiotest(cossim6870, num1, num2);
    
    % Feature: cosine similarity of nodes 49 and 69
    cossim4969 = cosine_similarity(:,49,69);
    [conmat_cossim4969(i,:,:), accuracy_cossim4969(i)] = likelihoodratiotest(cossim4969, num1, num2);
    
    % Feature: cosine similarity of nodes 66 and 68
    cossim6668 = cosine_similarity(:,66,68);
    [conmat_cossim6668(i,:,:), accuracy_cossim6668(i)] = likelihoodratiotest(cossim6668, num1, num2); 
    
    % Feature: average cosine similarity of all three pairs
    cossimall = mean(cosine_similarity(:,[68, 49, 66],[70, 69,68]),2);
    [conmat_cossimall(i,:,:), accuracy_cossimall(i)] = likelihoodratiotest(cossimall, num1, num2);  

end

mean_accuracy_deg6 = mean(accuracy_deg6);
mean_accuracy_deg20 = mean(accuracy_deg20);
mean_accuracy_deg44 = mean(accuracy_deg44);
mean_accuracy_degall = mean(accuracy_degall);
mean_accuracy_wdeg25 = mean(accuracy_wdeg25);
mean_accuracy_wdeg39 = mean(accuracy_wdeg39);
mean_accuracy_wdeg60 = mean(accuracy_wdeg60);
mean_accuracy_wdegall = mean(accuracy_wdegall);
mean_accuracy_eigcen23 = mean(accuracy_eigcen23);
mean_accuracy_eigcen25 = mean(accuracy_eigcen25);
mean_accuracy_eigcen36 = mean(accuracy_eigcen36);
mean_accuracy_eigcenall = mean(accuracy_eigcenall);
mean_accuracy_cossim6870 = mean(accuracy_cossim6870);
mean_accuracy_cossim4969 = mean(accuracy_cossim4969);
mean_accuracy_cossim6668= mean(accuracy_cossim6668);
mean_accuracy_cossimall = mean(accuracy_cossimall);

%% Plotting
accuracy_all = [mean_accuracy_deg6 mean_accuracy_deg20 mean_accuracy_deg44 mean_accuracy_degall; mean_accuracy_wdeg25...
mean_accuracy_wdeg39 mean_accuracy_wdeg60 mean_accuracy_wdegall; mean_accuracy_eigcen23 mean_accuracy_eigcen25 ...
mean_accuracy_eigcen36 mean_accuracy_eigcenall; mean_accuracy_cossim6870 mean_accuracy_cossim4969 ...
mean_accuracy_cossim6668 mean_accuracy_cossimall]*100;
figure(1)
h = bar(accuracy_all);

toplabels = {'#6';'#20';'#44';'Mean';'#25'; '#39'; '#60';'Mean';'#23';'#25';'#36';'Mean';'#68-70';'#49-69';'#66-68';'Mean'};
set(gca, 'XTickLabel',{'Degree','Weighted Degree','Eigenvector Centrality','Cosine Similarity'})
ybuff = 2;
yb = cat(1,h.YData);
xb = bsxfun(@plus,h(1).XData,[h.XOffset]');
hold on
for i = 1:size(yb,2)
    for j = 1:length(yb(:,1))    
        t = toplabels{4*i+j-4};
        text(xb(j,i),yb(j,i)+3,t,'VerticalAlignment','top','HorizontalAlignment','center')
    end
end
ylim([0 100])
ylabel('Accuracy [%]')
title('Accuracy of the Likelihood Ratio Classifier')
