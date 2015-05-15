clc; clear all; close all;

%%
%load data

%Get code from folder or dropbox(borrowed from gregs code/project tempalte)
fname = 'rohitg_statconn.zip';
files = dir();

clear found
for i=1:length(files)
    if strcmp(fname, files(i).name)
        found = 1;
        break;
    end
end

if ~exist('found', 'var')
    addr = ' https://www.dropbox.com/sh/idt3d0gylplyo31/AAAVUFyK3h6ZvSUTcQlLPyRaa/rohitg_statconn.zip';
    opts = ' --no-check-certificate';
    cmd = strcat('wget',' -O rohitg_statconn.zip',addr, opts);
    system(cmd, '-echo');
end

unzip(fname);

%~~~uncomment this part and comment the next part to run the second dataset~~~~
% load('mouse_cortex2.mat');
% A=mouse_cortex2;
% cell1=127;
% cell2=194-127;
% cell1_loc=[1:18 20:28 37:47 51:60 73:82 92:102 106:109 116:129 133:134 141:149 151:152 154:158 163:176 179:180 183:184 191:194]; 

load('mouse_cortex1.mat');
A=mouse_cortex1;
cell1=15; %number of the first type of cell
cell2=14; %number of the second type of cell
cell1_loc=[2:14 24:25]; %locations of the first type of cell

%%
%ordering locations of first cell type and then the second cell type
cell2_loc=[];
for i=1:(cell1+cell2)
    if isempty(find(cell1_loc==i))
        cell2_loc=[cell2_loc i];
    end
end
nodes=[cell1_loc cell2_loc];


%%
%Make a matrix, where edge indicated if there's a connection between
%neurons and analyze it. The matrix is separated into blocks.
M=zeros((cell1+cell2));
for i=1:length(A)
    M(nodes==A(i,1),nodes==A(i,2))=1;   
    M(nodes==A(i,2),nodes==A(i,1))=1;  
end

figure(1);hold all;
x=linspace(1,(cell1+cell2));y=cell1;
imagesc(M);plot(x,y,'k');plot(y,x,'k');colormap('winter');
axis([1 (cell1+cell2) 1 (cell1+cell2)]);
SBM=M;

%find estimated probability of cell type 2 connecting to cell type 1 that
%are connected/unconnected
[ P_est_2toc1,P_est_2tou1,P_est_2to1,P_est_1to1] = analyzemat( SBM,cell1,cell2 )

%run a permutation test to see the liklihood of the probability estimates
P_test=Ptest(cell1,cell2, P_est_2toc1,P_est_2tou1,P_est_2to1,P_est_1to1)