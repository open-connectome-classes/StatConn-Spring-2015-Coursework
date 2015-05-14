%Final Project - Class: Statistical Connectomics
%Author: Sandra Gomez R., May 2015
%Project: Comparison of SBM and biSBM for the infered C. elegans glia network

%% Should be submitted along with other files:

%Celegans_data.xlsx - the excel file with the data
%adjlist2matrix.m  - .m file
%.m files created by the Christopher Aicher group:
%wsbm.m - wrapper function for the weighted SBM
%main_alg.m - infers the data structure and feeds it to wsbm_driver
%wsbm_driver.m - outputs the best sbm for the data
%biwsbm.m - runs wsbm.m for bipartite networks

% Names of matrices used %
%GN - Glia-Neuron adj list
%NN - Neuron-Neuron adj list
%GG = Glia-Glia adj list

%% -------------------------------CODE----------------------------------- %%

%Load files - Celegans Xcel file
fprintf('loading data');
filename='Celegans_data.xlsx';
GN=xlsread(filename,2,'D4:E118');
NN=xlsread(filename,1,'G2:H6419');
GN_w=xlsread(filename,2,'D4:F118');
NN_w=xlsread(filename,1,'G2:I6419');%"weighted" versions; I'm cheating, the weighhts are all one
fprintf('Assigning variable names')
types_GN=xlsread(filename,5,'C1:C333');%these are the labels of the nodes; they're used in biSBM

%Download files? << STILL NOT WORKING ugh
%{
addr = 'http://sandragomez21.weebly.com/zip-files.html';
    opts = ' --no-check-certificate';
    system('C:\cygwin\bin\bash --login -c "addr"');
%}

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

Total=adjlist2matrix(Total);
figure (4); imagesc(Total+Total'); colorbar %Worm Glia+Neuron Connectome?
title('Worm Glia+Neuron Connectome');
figure(5); spy(Total);

%Get the "weighted" versions of the adjlist< wsbm needs a weighted adjlist
%third column is only ones...yeah, I'm cheating
%Total_w
total=ones(size(Total,1),1);
Total_w =[Total, total];
%GG_w
gg=ones(size(GG,1),1);
GG_w=[GG,gg];

%Run SBM for NN, GG and GN separately
fprintf('Running SBMs for networks');
%wsbm is running with a default k=4 clusters,
%that should be optimized, so I'll modify that later
fprintf('Running SBM for NN');
LabelsNN = wsbm(NN_w);
fprintf('Running SBM for GG');
LabelsGG = wsbm(GG_w);
fprintf('Running SBM for GN');
LabelsGN = wsbm(GN_w);
fprintf('Running SBM for Total');
LabelsTotal = wsbm(Total_w);

%Run bipartite SBM for GN_w

a=1; b=1; %these are the two k for the SBM..
[Labels1, Model1]=biwsbm(GN_w,a,b,types_GN);
%plot the Prob matrix
model1=Model1.Para.mu; %use only numeric part of statistical model
figure (6); imagesc(model1'); %Probability matrix
title('Probability of a vertex belonging to a group  GN_full matrix');
xlabel('Probability of belonging to group'); ylabel('Vertex');
colorbar

%Run bipartite SBM for Total_w

a=1; b=1; %these are the two k for the SBM..
[Labels2, Model2]=biwsbm(Total_w,a,b,types_total);
%plot the Prob matrix
model2=Model2.Para.mu; %use only numeric part of statistical model
figure (7); imagesc(model2'); %Probability matrix
title('Probability of a vertex belonging to a group - G+N WormConnectome');
xlabel('Probability of belonging to group'); ylabel('Vertex');
colorbar

