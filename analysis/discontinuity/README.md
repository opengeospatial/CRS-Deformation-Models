# Code to calculate okada fault model on a grid of points

Disclaimer: Everything here is untested.

## faultgrid.py

Script for calculating an Okada fault dislocation model on a grid.  This was set up for testing interpolation across a fault 
so sets up a grid, and then calculates on a subdivision of the grid both using Okada formulae and interpolation to see how 
they compare.

Input files are a fault definition file, see calc_okada.help for definition, and a grid definition file, see faultgrid.yaml for example.

Output files are based on name of the fault model definition file. 

calc_okada is an C++ program as defined in src.  Ubuntu (20.04) and windows .exe are included - they may work.

