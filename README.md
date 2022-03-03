# PSD

The code for calculating the value of the PSD verification measure

--- Description:

Calculates the value of the Precipitation Smoothing Distance (PSD) verification measure. PSD is a verification measure that tries to estimate the spatial displacement of precipitation in two fields. 

--- Usage:

function calculate_PSD(fa, fb)

--- Function arguments:

fa - a two-dimensional numpy array representing the first field. Only positive and zero values are allowed.

fb - a two-dimensional numpy array representing the second field. Only positive and zero values are allowed. The array's dimensions need to be the same as the dimensions of fa. 

--- Details:

The code utilizes the Fast-Fourier-Transform-based convolution for computationally efficient smoothing, along with the Bisection method. 

--- Return value:

The function returns a single numeric value representing the spatial distance (expressed as a number of grid points - see the example). 

In case of problems (e.g., if the two input arrays do not have the same dimensions, if the arrays contain negative values, ...), the function displays a warning message and returns value -1.

--- Required python libraries:

The code requires the following libraries to be installed: numpy, scipy, skimage

--- Example:

The example is in the file PSD_example.py.

--- Author:

Gregor Skok (Gregor.Skok@fmf.uni-lj.si)

--- References:

Paper submitted to Applied sciences journal. 

