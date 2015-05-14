# README

## About

This is Max C.'s StatConn final project. When executed, it will generate three figures in the Figures directory, showing ASE-MSE clustering's performance on

- An SBM
- The BA preferential attachment model
- BA applied to the C. elegans connectome

## Dependencies

You will need the following libraries (my version in parentheses). I recommend using pip (you will need a C++ compiler), and I've included the appropriate commands.

- Python 3 (3.4.3)
- NumPy (1.9.2)

	pip install numpy

- SciPy ( ... )

	pip install scipy

- Matplotlib ( ... )

	pip install matplotlib

- Scikit Learn ( ... )

	pip install sklearn

## Running

To run, simply execute

	python main.py

(or python3 if you have multiple installations). The relevant figures will appear in the Figures directory. (The random numbers at the end of the figure names ensure that no iinteresting old results get overwritten accidentally.)