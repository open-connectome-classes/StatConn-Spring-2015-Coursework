# README

## About

This is Max C.'s StatConn final project. When executed, it will generate three figures in the Figures directory, showing ASE-MSE clustering's performance on

- An SBM
- The BA preferential attachment model
- BA applied to the C. elegans connectome

## Dependencies

You will need the following standard numerical Python stack (my version in parentheses):

- Python 3     ( 3.4.1 )
- NumPy        ( 1.9.2 )
- SciPy        ( 0.14.1 )
- matplotlib   ( 1.4.3 )
- scikit-learn ( 0.15.2 )

## Running

To run, simply execute

	python main.py

(or python3 if you have multiple installations). The relevant figures will appear in the Figures directory. (The random numbers at the end of the figure names ensure that no interesting old results get overwritten accidentally.)

This will take *a while*, so please be patient. :)