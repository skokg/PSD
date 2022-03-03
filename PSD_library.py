import numpy as np
import scipy
import skimage
from skimage.draw import disk
from scipy import signal

# function that returns a circular kernel field
def construct_circular_kernel_field(r):
	r_int_max=np.int(np.floor(r))
	shape = (r_int_max*2+1, r_int_max*2+1)
	fcirc = np.zeros(shape)
	rr, cc = skimage.draw.disk((r_int_max, r_int_max), r+0.0001, shape=shape)
	fcirc[rr, cc] = 1
	return (fcirc/np.sum(fcirc))

def calculate_PSS(fa, fb, r, Q):
	
	fkernel=construct_circular_kernel_field(r=r)
	
	# smoothing - this function also automatically enlarges the domain area
	fa_smooth=scipy.signal.fftconvolve(fa, fkernel, mode='full')
	fb_smooth=scipy.signal.fftconvolve(fb, fkernel, mode='full')
	
	PSS=1.0-1.0/(2.0*float(fa.size)*Q)*(np.abs(fa_smooth - fb_smooth)).sum()
	
	return(PSS)

def calculate_PSD(fa, fb):
	
	# make a copy of the fields so the values in the input fields are not changed
	fa1=fa.copy()
	fb1=fb.copy()
	
	# compare dimensions of fields
	if fa1.shape != fb1.shape:
		print("ERROR: fa and fb arrays do not have the same shape. Returning -1 as result!")
		return(-1)
	
	# compare number of dimensions of fields
	if fa1.ndim != 2 or fb1.ndim != 2:
		print("ERROR: fa and fb arrays are not two-dimensional. Returning -1 as result!")
		return(-1)
	
	
	# compare the array has some elements
	if fa1.size == 0 or fb1.size == 0:
		print("ERROR: the dimensions of fa and fb arrays are zero. Returning -1 as result!")
		return(-1)
	
	# detect non-numeric values
	result=np.where(np.isfinite(fa1) == False)
	if len(result[0]) > 0:
		print("ERROR: fa or fb arrays contain some non-numeric values. Returning -1 as result!")
		return(-1)
	result=np.where(np.isfinite(fb1) == False)
	if len(result[0]) > 0:
		print("ERROR: fa or fb arrays contain some non-numeric values. Returning -1 as result!")
		return(-1)

	# detect masked array
	if isinstance(fa1, np.ma.MaskedArray) or isinstance(fb1, np.ma.MaskedArray) :
		print("ERROR: fa or fb arrays are masked arrays which is not allowed. Returning -1 as result!")
		return(-1)
	
	# detect negative values
	if fa1[fa1<0].size > 0 or fb1[fb1<0].size > 0 :
		print("ERROR: fa or fb arrays contain some negative values which is not allowed . Returning -1 as result!")
		return(-1)
	
	# cast to float - just in case the input fields are integers
	fa1=fa1.astype(float)
	fb1=fb1.astype(float)
	
	# Normalization
	fa1_avg=np.average(fa1)
	fb1_avg=np.average(fb1)
	
	if (fa1_avg == 0 or fb1_avg == 0):
		print("WARNING: At least one of the fields is empty (contains only zero values). Returning -1 as result!")
		return(-1)
	
	fa2=fa1.copy()
	fb2=fb1.copy()
	fa2=fa1/fa1_avg
	fb2=fb1/fb1_avg
	
	#special case if the two fields are indetical - return 0 in this case
	if (fa2==fb2).all():
		return(0)
	
	# Removal of overlapping values
	fa3=fa2.copy()
	fb3=fb2.copy()
	fa3=fa2 - np.minimum(fa2,fb2)
	fb3=fb2 - np.minimum(fa2,fb2)
	
	Q = fa3.sum()/fa2.sum()
	
	r1=1
	PSS1=calculate_PSS(fa=fa3, fb=fb3, r=r1, Q=Q)
	
	# special case when PSS > 0.5 at n=1
	if PSS1 > 0.5:
		return(1);
	
	# increase r in steps of 5 % 
	diagonal = np.sqrt(fa1.shape[0]*fa1.shape[0] + fa1.shape[1]*fa1.shape[1])
	dr = np.ceil(diagonal*0.05)
	
	r2=r1+dr
	PSS2=calculate_PSS(fa=fa3, fb=fb3, r=r2, Q=Q)
	while PSS2 < 0.5:
		r1=r2
		PSS1=PSS2
		r2=r1+dr
		PSS2=calculate_PSS(fa=fa3, fb=fb3, r=r2, Q=Q)
	
	# Bisection to get to the final r1 and r2 that bound the PSS=0.5 value
	while r2 - r1 > 1:
		# select middle point
		rnew=int((r1 + r2)/2);
		PSSnew=calculate_PSS(fa=fa3, fb=fb3, r=rnew, Q=Q)
		
		if PSSnew > 0.5:
			r2=rnew;
			PSS2=PSSnew;
		else:
			r1=rnew;
			PSS1=PSSnew;
	
	PSD = 0.808*Q*float(r2)
	
	return(PSD)

