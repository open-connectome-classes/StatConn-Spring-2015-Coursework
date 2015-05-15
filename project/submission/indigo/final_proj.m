%% "Clustering Algorithm Reverse Engineering of the Mouse Neocortex Connectome"
% Statistical Connectomics Final Project
% Written by: Indigo Vaux Loring Rose
% Designed in: MATLAB R2014b
% Due: Thursday May 14, 2015 at 23:59:59 EDT

clear; close all; clc;

%% Download Data
%{
%Download & unzip data using wget, as specified by Greg
system('wget https://www.dropbox.com/sh/idt3d0gylplyo31/AADfoAAgcWg-vq0SsRpx1u_la/Indigo.zip?dl=1');
system('unzip Indigo.zip*');

%Remove extra files from unzipping
system('rm Indigo.zip?dl=1');
system('rmdir __MACOSX');
%}
%% Import Data
ant_raw = csvread('out-ant.csv',1,1,[1,1,49,49]); % ant = anteriorgrade tracer data
ret_raw = csvread('out-ret.csv',1,1,[1,1,49,49]); % ret = retrograde tracer data

ant = weight_conversion(ant_raw, 'autofix'); %Makes data look cleaner to view
ret = weight_conversion(ret_raw, 'autofix');

%Read labels
formspec = '%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
lab_in = textscan(fopen('out-ant.csv','r'), formspec, 'Delimiter', ...
    ',',  'ReturnOnError', false);
lab_withPHAL = [lab_in{:, 1}];
lab = lab_withPHAL(2:end); %Excludes title line "PHAL"

%Import groundtruth data
%gt = csvread('Groundtruth.csv',1,1,[1,1,49,3]);

%% Plot Unclustered Data
figure(1); %subplot(1,2,1);
imagesc(ant); colorbar('FontSize',11);
title('Unclustered Anterograde Transport','FontSize',16);
ylabel('Cortical Area, Connecting From:','FontSize',14);
xlabel('Cortical Area, Connecting To:','FontSize',14);
set(gca, 'XTick', 1:49, 'XTickLabel', lab,'FontSize',8);
set(gca, 'YTick', 1:49, 'YTickLabel', lab,'FontSize',8);
rotateXLabels(gca(),90);

figure(2); %subplot(1,2,2);
imagesc(ret); colorbar('FontSize',11);                    
title('Unclustered Retrograde Transport','FontSize',16);
ylabel('Cortical Area, Connecting From:','FontSize',14);
xlabel('Cortical Area, Connecting To:','FontSize',14);
set(gca, 'XTick', 1:49, 'XTickLabel', lab,'FontSize',8);
set(gca, 'YTick', 1:49, 'YTickLabel', lab,'FontSize',8);
rotateXLabels(gca(),90);

%Using spy
%spy(ant);
%spy(ret);

%% Produce Composite Matrix

% Ideally, the weight of anteriograde transport from region A to region B
% should be the same as retrograde transport from region B to region A,
% since they should be using the same axon tracts

% Complile matricies into cmat. This has the same orientation as ant, so
% the regions listed in the rows map onto regions listed in the columns

cmat = ant + ret'; %unthresholded

% Threshold: Axon tracts that don't match back between ant and ret are set 
% to 0. This was specified in the extra material portion of the paper.

threshold = 2; %Anything <2 gets set to 0
cmat_thresh = threshold_absolute(cmat,threshold);

cmat_bin = weight_conversion(cmat_thresh, 'binarize'); 
%Talk about the implications of doing this in report

% Reverse engineer clustering algorithm used in paper
%% kmeans Algorithm 

subnet4 = 4; %The paper identified 4 major networks
subnet12 = 12; %It also identified 12 minor subnetwork model
clust_kmeans4 = kmeans(cmat_bin, subnet4);
clust_kmeans12= kmeans(cmat_bin, subnet12);

index = 1:49;
lab12 = [lab num2cell(index') num2cell(clust_kmeans12)];
labsort12 = sortcell(lab12,3);
lab4 = [lab num2cell(index') num2cell(clust_kmeans4)];
labsort4 = sortcell(lab4,3);

figure(3);
%plot(labsort4(2),labsort4(3));
subplot(2,1,1);
bar(cell2mat(labsort12(:,2)),cell2mat(labsort12(:,3)));
title('kmeans, k=12','FontSize',16);
xlabel('Brain Region Number','FontSize',14);

subplot(2,1,2);
bar(cell2mat(labsort4(:,2)),cell2mat(labsort4(:,3)));
title('kmeans, k=4','FontSize',16);
xlabel('Brain Region Number','FontSize',14);

%labn = cell2mat(lab);
%test = 1:49;
%{
A=1:49;
B=subnet4

for ii = 1:lenth(A)
     index = find(B == ii);
     C(ii) = A(index);
end

A = 1:49;
%A = lab;
B = clust_kmeans4';
C = [-B;A].';
sortout = sortrows(C);

sortout(:,2)
%}


figure(4);
subplot(1,2,1);
silhouette(cmat_bin,clust_kmeans4);
title('kmeans, k=4','FontSize',16);
ylabel('Cluster','FontSize',14);

%figure(4);
subplot(1,2,2);
silhouette(cmat_bin,clust_kmeans12);
title('kmeans, k=12','FontSize',16);
ylabel('Cluster','FontSize',14);


%% Heirarchical Clustering Algorithm 

clus = clusterdata(cmat_bin, 1);
dist = pdist(cmat_bin);
links = linkage(dist,'single');
linkc = linkage(dist,'complete');
linka = linkage(dist,'average');
linkw = linkage(dist,'weighted');

%figure; dendrogram(links);
figure(5); 
subplot(2,2,1);
dendrogram(links);
title('Dendrogram of Clusters Under Simple Linkage Method','FontSize',9);
xlabel('Neocortical Region Number');

subplot(2,2,2);
dendrogram(linkc);
title('Dendrogram of Clusters Under Complete Linkage Method','FontSize',9);
xlabel('Neocortical Region Number');
%set(gca, 'XTick', 1:49, 'XTickLabel', lab,'FontSize',8);
%rotateXLabels(gca(),90);

subplot(2,2,3);
dendrogram(linka); %These are being super weird
title('Dendrogram of Clusters Under Average Linkage Method','FontSize',9);
xlabel('Neocortical Region Number');

subplot(2,2,4);
dendrogram(linkw);
title('Dendrogram of Clusters Under Weighted Linkage Method','FontSize',9);
xlabel('Neocortical Region Number');


%% Rearrange Matrix

figure(6);
[order1 cmatsort]= reorder_mod(cmat,clust_kmeans12');
imagesc(cmatsort); colorbar('FontSize',11); 
xlabel('Brain Region Number','FontSize',14);
ylabel('Brain Region Number','FontSize',14);
set(gca, 'XTick', 1:49, 'XTickLabel', order1,'FontSize',8);
set(gca, 'YTick', 1:49, 'YTickLabel', order1,'FontSize',8);
title('Ordered Combined Matrix','FontSize',16);
%{
cmat = cmat_bin + cmat_bin';
r = symrcm(cmat);
p = symamd(cmat);
figure; imagesc(cmat);
figure; imagesc(cmat(r,r));
figure; imagesc(cmat(p,p));
%}



% Compare together
