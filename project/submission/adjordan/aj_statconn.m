% % % AJ_STATCONN.M
% % My final project for statconn compares modularity values of preictal and
% % ictal networks for seizure sEEG data.
% % Austin Jordan; May 14th, 2015
% % 
% % The script must be run from a directory structured as follows:
% % 
% % -Current Directory:
% %     -aj_statconn.m
% %     -cluster_jl.m
% %     -getAllFiles.m
% %     -Modularity_Cluster_Seizure.m
% % 
% % ... only external dependency is wget

% % % Preamble
% % Check if zip file has already been downloaded
fname = 'austin_jordan_data.zip';
files = dir();

clear found
for i = 1:length(files)
    if strcmp(fname, files(i).name)
        found = 1;
        break;
    end
end

% % If it hasn't already been downloaded...
if ~exist('found','var')
    addr = ' https://www.dropbox.com/sh/idt3d0gylplyo31/AACtO4G3GPq9dRu7a2ggzGHla/austin_jordan_data.zip';
    opts = ' --no-check-certificate';
    cmd = strcat('wget',addr,opts);
    system(cmd, '-echo');
end

% % Unzip the file
unzip(fname)

% % % Preamble over. Enter my code.
clear all
close all

% % Get a list of all of the paths to the .mat files I will be using
file_list = getAllFiles(pwd);

% % Count how many there are
num_files = size(file_list,1);

% % Look through each file to use with my code, Modularity_Cluster_Seizure
for i = 1:num_files
    file = file_list{i};
    if strcmp('/Users/AustinJordan/StatConn/Final_Project/EZT053_seiz01.mat',file)
        plot = 1; % For plotting representative figure
    end
    load(file)
    Modularity_Cluster_Seizure(preictal_adj_mat,ictal_adj_mat,plot);
    plot = 0;
end