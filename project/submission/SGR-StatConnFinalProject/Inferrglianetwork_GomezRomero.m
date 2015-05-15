%Final Project - Class: Statistical Connectomics

%Author: Sandra Gomez R., May 2015
%Software: Created on MATLAB R2014b
%Project: Clustering and inferring the C. elegans glia network

%% Description : %%

%This script uses the data from the WormAtlas project to inferr the C. elegans
%glial cell connectivity matrix: It takes the neuron-neuron adjacecy list
%and neuron-glia adjacency list and inferrs the glia-glia adjacency list
%under the assumption that glial cells that are connected to a given neuron A
%are connected to the glial cells connected to neuron B given that neuron A
%and B are connected. 

%The inferred glia-glia network is called GG, and added to the
%neuron-neuron adjacency list and the neuron-glia adjacency list to create
%the matrix 'Total', which we call the 'Worm Glia+Neuron Connectome'
%because it includes all connections.

%It then performs clustering on 'Total' and on 'GN' (the glia-neuron
%matrix) for comparison.

%% INPUT %%

%Celegans_data.xlsx - the excel file with the data; the script downloads
%it and unzips it

%% SUBFUNCTIONS %%

%adjlist2matrix.m  - .m file to convert form adjacency list to matrix

%% VARIABLES NAMES %%

% Names of adjacency lists used %
%GN - Glia-Neuron adj list
%NN - Neuron-Neuron adj list
%GG = inferred Glia-Glia adj list
%Total = Neuron-neuron, neuron-glia and glia-glia adj list

%% -------------------------------CODE----------------------------------- %%

% Download the data << copied from Greg's code

fname = 'SGR-StatConnFinalProject.zip';
files = dir();

clear found
for i=1:length(files) %silly but your directory should be almost empty
    if strcmp(fname, files(i).name)
        found = 1;
        break;
    end
end

%if no, download it
if ~exist('found', 'var')
    %craft and run command
    addr = ' https://www.dropbox.com/sh/idt3d0gylplyo31/AADFjehGbjn3hkmBsbl3Lmuaa/SGR-StatConnFinalProject.zip';
    opts = ' --no-check-certificate';
    cmd = strcat('wget', addr, opts);
    system(cmd, '-echo');
end

%Load files - Celegans Xcel file
fprintf('loading data');
filename='Celegans_data.xlsx';
GN=xlsread(filename,2,'D4:E118');
NN=xlsread(filename,1,'G2:H6419');
GN_w=xlsread(filename,2,'D4:F118');
NN_w=xlsread(filename,1,'G2:I6419');%"weighted" versions; I'm cheating, the weighhts are all one
fprintf('Assigning variable names')
types_GN=xlsread(filename,5,'C1:C333');%these are the labels of the nodes; they're used in biSBM

%Infer Glia-Glia adj list
fprintf('Constructing Glia-Glia Network');
T=unique(GN(:,2)); %unique neurons
n=size(GN,1);
GG=zeros(1,2);%lame preallocation
for i=1:n
    [~,locNN]=ismember(GN(i,2),NN(:,1));%get the position in NN of the match
    [~,loct]=ismember(NN(locNN,1),T(:,1));
    a=NN(NN(:,1)==T(loct,1),2);%neurons connected to T(loct,1)
    b=GN(T(loct,1)==GN(:,2),1);%glia connected to T(loct,1)
    l=size(a,1);
    store=zeros(1,1);
    for k=1:l
        
        t=GN(a(k)==GN(:,2),1);%glia connected to a(k)
        if isempty(t)
            t=[];
        else
            store=[store;t];
        end
        glia=[b;store]; %these are all the glia that touch T(loct,1) in some way
        glia(glia==0)=[];%get rid of those nasty zeroes
        glia=unique(glia); %avoid repetition
        lil_gg=combnk(glia,2); %it's a tiny glia-glia adj list
    end
    GG=[GG;lil_gg];
end
GG(1,:)=[];%bye zeroes << so in this adj list, glia are named from 282 to 331

fprintf('Obtaining Total G+N matrix');
%Mix the connection data from original GN, NN networkw with GG
Total=[GG ; GN; NN ];

%create edge type labels for future use
fprintf('Creating data labels for lists');
w=size(GG,1);f=size(GN,1); q=size(NN,1);
types_total=ones(size(Total,1),1);
types_total(1:w,1)=types_total(1:w,1)*1;
types_total(w+1:w+f,1)=types_total(w+1:w+f,1)*2;
types_total(f+1:f+q,1)=types_total(f+1:f+q,1)*3;

%Obtain the only-glia network, so rename glia form 1 to 50
GG=GG-281; %substract 281


%Obtain and visualize full matrices
fprintf('Converting adjacency lists to matrices');

NN_full=adjlist2matrix(NN);
figure (1); imagesc(NN_full+NN_full'); colorbar %Worm Neuron Connectome
title('Worm Connectome - Neurons');

GG_full=adjlist2matrix(GG);
figure (2); imagesc(GG_full+GG_full'); colorbar %Worm Glia Connectome?
title('Worm Glia Connectome');

GN_full=adjlist2matrix(GN);
figure (3); imagesc(GN_full+GN_full'); colorbar %Worm Glia+Neuron Connectome?
title('Worm GN matrix Graph');

Total_full=adjlist2matrix(Total);
figure (4); imagesc(Total_full+Total_full'); colorbar %Worm Glia+Neuron Connectome?
title('Worm Glia+Neuron Connectome');

figure(5); spy(Total_full);
title('Visualizing the Worm Glia+Neuron Connectome');

%Clustering analysis

%Kmeans 

fprintf('Running kmeans');
%For the 'Total' matrix; 
%this is pretty obvious. There are two clusters: glia and neurons
%Run k means for Total
Test1=kmeans(Total_full,2);
figure (6); silhouette(Total_full,Test1); %visaulize
title('keamns for Worm Glia+Neuron - 2 clusters');
%Run kmeans on GN for comparison 
Test2=kmeans(GN_full,2);
figure (7); silhouette(GN_full,Test2); %visualize
title('kmeans for Worm GN matrix Graph - 2 clusters');

%Hierarchical clustering

%For the 'Total' matrix
%we're gonna cluster for 2 and 3 clusters.
fprintf('hierarchical clustering of Total matrix');
Cluster1_total=clusterdata(Total,'maxclust',2);% 2 clusters
figure (8); scatter(Total(:,1),Total(:,2),10,Cluster1_total, 'filled');
title('Hierarchical Clustering of Glia+Nueron connectome -2 clusters');

Cluster1_total=clusterdata(Total,'maxclust',3);% 3 clusters
figure (9); scatter(Total(:,1),Total(:,2),10,Cluster1_total, 'filled');
title('Hierarchical Clustering of Glia+Nueron connectome -3 clusters');

%For the GN matrix for comparison
fprintf('hierarchical clustering of GN matrix');
Cluster1_gn=clusterdata(GN,'maxclust',2);%2 clusters
figure (10); scatter(GN(:,1),GN(:,2),10,Cluster1_gn, 'filled');
title('Hierarchical clustering of GN matrix - 2 clusters');

Cluster2_gn=clusterdata(GN,'maxclust',3);%3 clusters
figure (11); scatter(GN(:,1),GN(:,2),10,Cluster2_gn, 'filled');
title('Hierarchical clustering of GN matrix - 3 clusters');

%Analyzing GG - Now, let's study how the Astrocyte network is organized

%Analyzing GG
fprintf('Analyzing the structure of teh glia-glia network');
Test_gg1=kmeans(Total_full,2);
figure (12); silhouette(Total_full,Test_gg1); %visualize
title('kmeans for Worm glia matrix Graph - 2 clusters');

Cluster1_gg=clusterdata(GG,'maxclust',2);%2 clusters
figure (13); scatter(GG(:,1),GG(:,2),10,Cluster1_gg, 'filled');
title('Hierarchical clustering of Worm Glia connectome - 2 clusters');

Cluster2_gg=clusterdata(GG,'maxclust',3);%3 clusters
figure (14); scatter(GG(:,1),GG(:,2),10, Cluster2_gg, 'filled');
title('Hierarchical clustering of Worm Glia connectome - 3 clusters');



