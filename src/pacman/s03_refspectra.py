import numpy as np
import os
from scipy.interpolate import interp1d
from ..lib import plots
from ..lib import stellar_spectrum
from ..lib import manageevent as me
from ..lib import util


def run03(eventlabel, workdir, meta=None):
    """
    Retrieves the bandpass (G102 or G141) and the stellar spectrum and takes the product to create a reference spectrum.

    Options for the stellar model:
    - Blackbody
    - k93models
    - ck04models
    - phoenix

    The last three stellar models are retrieved from https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/

	Parameters
	----------
	eventlabel : str
	   the label given to the event in the run script. Will determine the name of the run directory
	workdir : str
	   the name of the work directory.
	meta
	   the name of the metadata file

	Returns
	-------
	meta
	   meta object with all the meta data stored in s02

    Notes:
    ----------
    History:
        Written by Sebastian Zieba      December 2021
    """

    print('Starting s03')

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")


    ### Stellar Spectrum
    Teff, logg, MH = meta.Teff, meta.logg, meta.MH
    print('Using {0} model.\n'.format(meta.sm))

    if meta.sm in ['k93models', 'ck04models', 'phoenix']:
        sm_wvl, sm_flux = stellar_spectrum.get_sm(meta, MH, logg, Teff)
    elif meta.sm == 'blackbody':
        sm_wvl, sm_flux = stellar_spectrum.get_bb(Teff)
    else:
        print('You have not entered a valid stellar model.\n'
              'Options are k93models, ck04models, phoenix or blackbody \n'
              'We proceed with the bb')
        sm_wvl, sm_flux = stellar_spectrum.get_bb(Teff)

    #only store the spectrum between 0.1 microns and 10 microns
    sm_wvl_mask = np.bitwise_and(0.1e-6 < sm_wvl, sm_wvl < 10e-6)
    sm_wvl = sm_wvl[sm_wvl_mask]
    sm_flux = sm_flux[sm_wvl_mask]

    if meta.smooth:
        sm_wvl, sm_flux = util.gaussian_kernel(meta, sm_wvl, sm_flux)


    ### Bandpass
    if meta.grism == 'G141':
        grism = 'g141'
    elif meta.grism == 'G102':
        grism = 'g102'
    print('Using {0} grism.'.format(grism))

    #Read in bandpass for the used grism
    bp_wvl, bp_val = np.loadtxt(meta.pacmandir + '/ancil/bandpass/bandpass_{0}.txt'.format(grism)).T


    ### Creating the reference spectrum
    bp_wvl = bp_wvl * 1e-10
    bp_val = bp_val / max(bp_val)
    #sm_flux = sm_flux #* sm_wvl  # in order to convert from W/m^3/sr units to W/m^2/sr
    #sm_flux = sm_flux / max(sm_flux)

    meta.refspecdir = meta.workdir + '/ancil/refspec/'
    if not os.path.exists(meta.refspecdir):
        os.mkdir(meta.refspecdir)

    #Interpolate stellar model so that we can multiply it with the bandpass
    f = interp1d(sm_wvl, sm_flux, kind='linear')

    ref_wvl = bp_wvl
    ref_flux = f(bp_wvl) * bp_val
    ref_flux = ref_flux / max(ref_flux)

    # Save reference spectrum
    np.savetxt(meta.refspecdir + '/refspec.txt', list(zip(ref_wvl, ref_flux)))

    if meta.save_refspec_plot or meta.show_refspec_plot:
        plots.refspec(bp_wvl, bp_val, sm_wvl, sm_flux, ref_wvl, ref_flux, meta)

    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    print('Finished s03 \n')

    return meta
