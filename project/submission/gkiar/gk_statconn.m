%% GK_STATCONN.M
% My final project for statconn compares TRT performance of different
% parcelation schemes for human brain graphs.
% Greg Kiar; Apr 29th, 2015
%
%The script must be run from in a directory structured as follows:
% **(w/ nothing else)**
% *(The reason for this is so I can cleanly read in the 168 graphs...)*
%
%-Current Directory:
%   -gk_statconn.m
%   -eval_partitions.m
%   [the following are contents of a ZIP this script will download]
%   -gk_statconn.zip
%   -subject_information.csv
%   -partition1/:
%       -graph1.mat
%       -...
%       -graphN.mat
%   -....
%   -partitionN/:
%       -graph1.mat
%       -...
%       -graphN.mat
%
% ... only external dependency is wget

%% Preamble

%have I downloaded my zip file
fname = 'gk_statconn.zip';
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
    addr = ' https://www.dropbox.com/sh/idt3d0gylplyo31/AADFjehGbjn3hkmBsbl3Lmuaa/gk_statconn.zip';
    opts = ' --no-check-certificate';
    cmd = strcat('wget', addr, opts);
    system(cmd, '-echo');
end

%regardless, unzip it
unzip(fname);

%% Evaluate
clear performance
performance = eval_partitions();

%% Plot?
% currently happening inside eval_partitions()...