
Final Project - Class: Statistical Connectomics

Author: Sandra Gomez R., May 2015
Software: Created on MATLAB R2014b
Project: Inferring and clustering the C. elegans glia network

%% To run this file %%
Simply run the 'Inferrglianetwork_GomezRomero.m' file. It will download the Excel file
'Celegans_data.xlsx' and call the subfunction 'adjlist2matrix'.

%% Description %%

This script uses the data from the WormAtlas (hermaphrodite neuron-neuron adjacency 
list and neuron-glia adjacency list) project to inferr the hermaphrodite C. elegans
glial cell connectivity matrix. It works under the assumption that glial cells 
connected to neuron A are connected to glial cells connected to neuron B given 
that neuron A and B are connected. 

The inferred glia-glia network is then called GG, and added to the
neuron-neuron  (NN) adjacency list and the neuron-glia (GN) adjacency list to create
the matrix named 'Total', which we call the 'Worm Glia+Neuron Connectome'
because it includes all types of connections between the two cell types.

The script then performs kmeans and hierarchical clustering on 'Total' and on 'GN' 
(the glia-neuron matrix) for comparison. 

%% INPUT %%

Celegans_data.xlsx - the excel file with the data; the script downloads it and unzips it

%% Outputs %%

GG = Glia-glia adj list
GG_full = Glia-glia adj matrix
Total = glia-neuron-neuron adj list
Total_full= glia-neuron-neuron adj matrix

Lost of figures. 14 in total

%% SUBFUNCTIONS %%

adjlist2matrix.m  - .m file to convert form adjacency list to matrix

%% VARIABLES NAMES %%

% Names of adjacency lists used %

GN - Glia-Neuron adj list
NN - Neuron-Neuron adj list
GG = inferred Glia-Glia adj list
Total = Neuron-neuron, neuron-glia and glia-glia adj list