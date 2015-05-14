% Alan Juliano
% Homework 1


% The goal of this assignment is to create an example when the kmeans task
% fails. The Kmeans process is a clustering techniqe relating to vector
% quantization. This is very important in data mining when doing
% classification techniques. The n observations are mined into clusters
% relying on the nearest mean as the prototype of the cluster, thus
% creating partitions within the data space. 

% The following data set would typically be interpreted by implementing
% seperate clusters for A/B and D/E/F. C, interestingly, could be placed in
% either cluster. K means, conversely, can not be implemented on this data.


    
%    A B C D E F     
a = [0,1,1,0,0,0; %A
     1,0,1,0,0,0; %B
     1,1,0,1,1,1; %C
     0,0,1,0,1,1; %D
     0,0,0,1,0,1; %E
     0,0,0,1,1,0]; %F
 
 
 
imagesc(a)

for i = 1:6
    kmeans(a,2)
end
 
%  This code is does not implement kmeans in an effective way
% and fails. Specifically, C/E/F nodes are one cluster while the other nodes
% are all clustered together. K means fails to implement itself on this
% data. 

