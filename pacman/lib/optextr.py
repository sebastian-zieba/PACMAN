import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import numpy.ma as ma


#Applies a boxcar smooth of length nsmooth to the vector x
#returns the smoothed vector
def smooth(x, nsmooth):
	#interpolates over masked values
	if (sum(x==0)> 0)&(sum(x)>0):
		bpix = x==0.0
		gpix = ~bpix
		gx = x[gpix]
		interp = np.interp(bpix.nonzero()[0], gpix.nonzero()[0], gx)
		x[bpix] = np.float32(interp)				
	return scipy.signal.medfilt(x, nsmooth)				#median filters the data

def diagnostics_plot(D, M, indmax, outlier_array, f_opt, profile):
	indmax = np.argmax(outlier_array)			#finds biggest outlier
	indmax = unravel_index(indmax, outlier_array.shape)	#converts it from flat to tuple

	plt.subplot(221)
	plt.title("Raw Data")
	plt.imshow(D,vmin=0, vmax=50)
	plt.scatter(x = indmax[1], y = indmax[0], color='w', marker='x')
	m = cm.ScalarMappable(cmap=  cm.jet)
	m.set_array(D)
	plt.colorbar(m)
	plt.subplot(222)
	plt.title("Outliers")
	plt.imshow(M*outlier_array, vmin=0, vmax=20)
	plt.scatter(x = indmax[1], y = indmax[0], color='w', marker='x')
	m.set_array(M*outlier_array)
	plt.colorbar(m)
	plt.subplot(222)
	plt.subplot(223)
	plt.title("Cut in spatial direction")
	plt.axvline(x = indmax[0], color="red")
	plt.plot(D[:, indmax[1]], label = "data")
	plt.plot((f_opt*profile)[:, indmax[1]], color="orange", label= "model")
	plt.legend()
	plt.xlabel('Wavelength [um]')
	plt.ylabel('Counts')
	plt.subplot(224)
	plt.title("Outliers: cut in spatial direction")
	plt.plot(outlier_array[:,indmax[1]])
	plt.axvline(x = indmax[0], color="red")
	plt.ylabel('Residuals')
	plt.xlabel('Wavelength [um]')
	plt.tight_layout()

	plt.show()
	plt.clf()

"""Function to optimally extract a spectrum:
Inputs:
	D: data array (already background subtracted)
	err: error array (in addition to photon noise; e.g. error due to background subtraction)
	f_std: box-extracted spectrum (from step 4 of Horne)
	var_std: variance of standard spectrum (also from step 4)
	M: array masking bad pixels; 0 is bad and 1 is good
	nsmooth: number of pixels to smooth over to estimate the spatial profile (7 works well)
	sig_cut: cutoff sigma for flagging outliers (10.0 works well)
	diagnostics: boolean flag specifying whether to make diagnostic plots
outputs:
	f_opt, var_opt:  optimally extracted spectrum and its variance"""

def optextr(D, err, f_std, var_std, M, nsmooth, sig_cut, diagnostics):
	#STEPS 5-8:  estimating spatial profile and removing cosmic rays
	f_opt = np.copy(f_std)							#array to store the optimally extracted spectrum
	outlier_array = np.zeros_like(D)					#array used to find outliers
	outliers = True 								
	numoutliers = 0								#number of outliers rejected by optimal extraction
	while outliers == True:
		#STEP 5:  construct spatial profile
		#interpolate over masked regions to better estimate spatial profile
		profile = np.apply_along_axis(smooth, 1, M*D, nsmooth)	

		#enforce positivity 
		ind = profile < 0.0
		profile[ind] = 0.0  			

		#handles case where whole column is 0
		ind = np.where(profile.sum(axis=0)==0)
		profile[:, ind] = 1.0				#since we normalize by column, we just need to set all the whole column equal to the same (nonzero) value 

		#normalize
		profile = profile/profile.sum(axis = 0)
    
		
		#STEP 6:  revise variance estimates
		var = abs(f_opt*profile) + err 

		#STEP 7:  mask cosmic rays/bad pixels
		outlier_array = M*(D - (f_opt*profile))**2/var		#number of standard deviations away from expected is each pixel

		maxes = np.argmax(outlier_array, axis = 0)
		#print maxes
		newoutliers = 0
		for ii in range(len(maxes)):
			ind2, ind1 = ii, maxes[ii]
			#print ind1, ind2
			if outlier_array[ind1, ind2] > sig_cut**2.:
				M[ind1, ind2] = 0.0
				numoutliers += 1
				newoutliers += 1
		if newoutliers == 0: outliers = False

		indmax = np.argmax(outlier_array)			#finds biggest outlier
		indmax = unravel_index(indmax, outlier_array.shape)	#converts it from flat to tuple

		#if outlier_array[indmax] > sig_cut**2.0: 
		#	M[indmax] = 0.0					#checks to see if the pixel is an outlier > sig_cut, and if so, masks that pixel
		#	numoutliers += 1
		#else: outliers = False					#if not outliers, switches outliers flag to false to close loop

		#STEP 8:  extract optimal spectrum
		f_opt = ((M*profile*D/var).sum(axis = 0))/(M*profile**2/var).sum(axis=0) 
	
	var_opt = (M*profile).sum(axis = 0)/(M*profile**2/var).sum(axis = 0)

	if diagnostics == True: diagnostics_plot(D, M, indmax, outlier_array, f_opt, profile)

	return f_opt, var_opt, numoutliers


