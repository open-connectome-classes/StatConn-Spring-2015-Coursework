This is the readme for the final Statistical Connectomics Project for 
William Gray Roncal.

This project explores the sources of graph error that result from the 
best available automated EM graph reconstructions, and possible methods
to reduce this error.

This submission consists of the following:

- this readme.txt file
- wrgr_statconn_driver.m:  Run this function to generate the overall results
- cleanup_speckle.m:  A utility function to cleanup source data
- statconn_compute_graph_error_simple.m:  A function that computes graph error,
    given a test and true line graph
- statconn_construct_graph.m:  Function to build neuron and line graphs from edge lists
- statconn_associate_neusyn.m: Estimates graph by associating neuron and synapse label data

This submission also depends on the CAJAL toolbox, which needs to be setup prior to running:

https://github.com/openconnectome/cajal
(to initialize the toolbox, follow the instructions, typically:
run('<toolbox root>/cajal3d.m')

The driver script downloads data from both the Open Connectome Project and Dropbox.

