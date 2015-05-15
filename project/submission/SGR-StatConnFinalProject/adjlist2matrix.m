%Final Project - Class: Statistical Connectomics
%Author: Sandra Gomez R., May 2015
%Software: Created on MATLAB R2014b
%Project: Clustering and inferring the C. elegans glia network

function [B]=adjlist2matrix(A)

%Returns the matrix of any adjacency list where the 
%input is A=any adj list with 2 columns 
rows = A(:,1);
cols = A(:,2);
s=ones(size(A,1),1);
BA = sparse([rows; cols],[cols; rows],[s,s]);
B=full(BA);
end