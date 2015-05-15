quick readme

code is in the code directory

downloading will use Python's urlib2 module and put the data into the code/data
directory.

Running
========
Please cd into the code directory before running.

Everything is written in Python 2.

$ python project.py

This will run the project, download data, and run parameter optimization to try
to fit a Random Dot Product Graph model to the C. Elegans connectome.

Required Python libraries
==================

networkx
scipy
scipy.optimize
scipy.stats
numpy
matplotlib
pylab (I think this is included in scipy
cma for CMA-ES optimizer

For parsing Excel data from the worm atlas (Not necessary -- you will be using
the derived file that I created from this -- but I included my parsing script 
anyway...)
xlrd

Code
=======

