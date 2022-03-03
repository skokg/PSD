# ---------------------------------------------
# A simple example with two 500 x 500 fields
# ---------------------------------------------

# import PSD functions from the library file
from PSD_library import *

#set domain size
dimx=500

# generate two fields filled with zeroes
fa=np.zeros((dimx, dimx))
fb=fa.copy()

# in the fields define two identicaly sized (20 x 20 grid points) non-zero regions with values=1 that are diagonaly displaced by 100 grid points. 
fa[200:220,200:220]=1
fb[200:220,300:320]=1

# calulate PSD value
PSD=calculate_PSD(fa, fb)

# print PSD value 
print(PSD)

# The example should output 100.1920. which means that the spatial displacement estimated by the PSD is 100.1920 grid points.


