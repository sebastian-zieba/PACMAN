import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal


def smooth(x, nsmooth):
    """Applies a boxcar smooth of length nsmooth to the vector x.
    Returns the smoothed vector."""
    #interpolates over masked values
    if (sum(x == 0) > 0) & (sum(x) > 0):
        bpix = x == 0.0
        gpix = ~bpix
        gx = x[gpix]
        interp = np.interp(bpix.nonzero()[0], gpix.nonzero()[0], gx)
        x[bpix] = np.float32(interp)
    return scipy.signal.medfilt(x, nsmooth)				#median filters the data


def diagnostics_plot(D, M, indmax, outlier_array, f_opt, profile, i, ii, meta):
    indmax = np.argmax(outlier_array)			#finds biggest outlier
    indmax = np.unravel_index(indmax, outlier_array.shape)	#converts it from flat to tuple
    fig,ax = plt.subplots(2, 2, figsize=(10, 10))

    ax[0,0].set_title("Raw Data")
    im00 = ax[0, 0].imshow(D, vmin=0, vmax=50)
    ax[0, 0].scatter(x=indmax[1], y=indmax[0], color='w', marker='x')
    cbar00 = fig.colorbar(im00, ax=ax[0, 0], shrink=0.89, pad=0.02)
    #m = plt.cm.ScalarMappable(cmap=  plt.cm.jet)
    #m.set_array(D)
    #ax[0,0].set_colorbar(m)

    ax[0,1].set_title("Outliers")
    im01 = ax[0,1].imshow(M*outlier_array,vmin=0, vmax=20)
    ax[0,1].scatter(x=indmax[1], y=indmax[0], color='w', marker='x')
    cbar01 = fig.colorbar(im01, ax=ax[0, 1], shrink=0.89, pad=0.02)
    #plt.colorbar(m)

    ax[1, 0].set_title("Cut in spatial direction")
    ax[1, 0].axvline(x=indmax[0], color="red")
    ax[1, 0].plot(D[:, indmax[1]], label="data")
    ax[1, 0].plot((f_opt*profile)[:, indmax[1]], color="orange", label="model")
    ax[1, 0].legend()
    ax[1, 0].set_xlabel('Wavelength [um]')
    ax[1, 0].set_ylabel('Counts')

    ax[1, 1].set_title("Outliers: cut in spatial direction")
    ax[1, 1].plot(outlier_array[:, indmax[1]])
    ax[1, 1].axvline(x=indmax[0], color="red")
    ax[1, 1].set_ylabel('Residuals')
    ax[1, 1].set_xlabel('Wavelength [um]')

    plt.tight_layout()
    opextr_dir = meta.workdir / 'figs' / 's20_optextr'
    if not opextr_dir.exists():
        opextr_dir.mkdir(parents=True)
    plt.savefig(opextr_dir / f'optextr{i}-{ii}.png', ii,
                dpi=120, bbox_inches='tight', pad_inches=0.05)
    plt.close('all')
    plt.clf()


def optextr(D, err, f_std, var_std, M, nsmooth,
            sig_cut, save_optextr_plot, i_sp, ii_sp, meta):
    """Function to optimally extract a spectrum.

    Parameters
    ----------
    D:
        data array (already background subtracted)
    err:
        error array (in addition to photon noise; e.g. error due to background subtraction)
    f_std:
        box-extracted spectrum (from step 4 of Horne)
    var_std:
        variance of standard spectrum (also from step 4)
    M:
        array masking bad pixels; 0 is bad and 1 is good
    nsmooth:
        number of pixels to smooth over to estimate the spatial profile (7 works well)
    sig_cut:
        cutoff sigma for flagging outliers (10.0 works well)
    diagnostics:
        boolean flag specifying whether to make diagnostic plots

    Returns
    -------
    f_opt, var_opt:
        optimally extracted spectrum and its variance
    """
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
        indmax = np.unravel_index(indmax, outlier_array.shape)	#converts it from flat to tuple

        #if outlier_array[indmax] > sig_cut**2.0:
        #	M[indmax] = 0.0					#checks to see if the pixel is an outlier > sig_cut, and if so, masks that pixel
        #	numoutliers += 1
        #else: outliers = False					#if not outliers, switches outliers flag to false to close loop

        #STEP 8:  extract optimal spectrum
        f_opt = ((M*profile*D/var).sum(axis = 0))/(M*profile**2/var).sum(axis=0)
        #(M*profile)[:,0] = 0 breaks everything!
        #((M*profile)).sum(axis = 0)=[0,1,1,1,1,...,1,1,1,1]
        #M[:,0]*profile[:,0] = [0,0,0,0,0]
        #       M[:,0] = [1,0,1,1,1,1,....,1,1,1,1]
        # profile[:,0] = [0,1,0,0,0,0,....,0,0,0,0]
        # problem occurs if eg there is only one non zero value in profile and M is zero where profile is non zero

    var_opt = (M*profile).sum(axis = 0)/(M*profile**2/var).sum(axis = 0)

    if save_optextr_plot == True: diagnostics_plot(D, M, indmax, outlier_array, f_opt, profile, i_sp, ii_sp, meta)

    return f_opt, var_opt, numoutliers
