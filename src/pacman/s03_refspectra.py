from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

from .lib import plots
from .lib import stellar_spectrum
from .lib import manageevent as me
from .lib import util


def run03(eventlabel, workdir: Path, meta=None):
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
        meta = me.loadevent(workdir / f'WFC3_{eventlabel}_Meta_Save')

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

    # UNITS: wavelength: nm

    #Keeping the whole spectral range of the spectrum is not neccessary
    # Let's only keep the spectrum around the G102 and G141 grisms, so let's say like 0.3 microns and 2.1 microns
    sm_wvl_mask = np.bitwise_and(300e-9 < sm_wvl, sm_wvl < 2.10e-6)
    sm_wvl = sm_wvl[sm_wvl_mask]
    sm_flux = sm_flux[sm_wvl_mask]

    # Let's smooth the stellar model using a gaussian kernel
    # More information here: https://ui.adsabs.harvard.edu/abs/2013ApJ...774...95D/abstract
    if meta.smooth:
        sm_wvl, sm_flux = util.gaussian_kernel(meta, sm_wvl, sm_flux)

    ### Bandpass
    grism = ""
    if meta.grism == 'G141':
        grism = 'g141'
    elif meta.grism == 'G102':
        grism = 'g102'
    print(f'Using {grism} grism.')

    #Read in bandpass for the used grism
    bp_wvl, bp_val = np.loadtxt(meta.pacmandir / 'data' / 'bandpass' / f'bandpass_{grism}.txt').T

    ### Creating the reference spectrum
    bp_wvl = bp_wvl * 1e-10
    bp_val = bp_val / max(bp_val)
    #sm_flux = sm_flux #* sm_wvl  # in order to convert from W/m^3/sr units to W/m^2/sr
    #sm_flux = sm_flux / max(sm_flux)

    meta.refspecdir = meta.workdir / 'ancil' / 'refspec'
    if not meta.refspecdir.exists():
        meta.refspecdir.mkdir(parents=True)

    #Interpolate stellar model so that we can multiply it with the bandpass
    f = interp1d(sm_wvl, sm_flux, kind='linear')

    ref_wvl = bp_wvl
    ref_flux = f(bp_wvl) * bp_val
    ref_flux = ref_flux / max(ref_flux)

    # Save reference spectrum
    np.savetxt(meta.refspecdir / 'refspec.txt', list(zip(ref_wvl, ref_flux)))

    if meta.save_refspec_plot or meta.show_refspec_plot:
        plots.refspec(bp_wvl, bp_val, sm_wvl, sm_flux, ref_wvl, ref_flux, meta)

    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir / f'WFC3_{meta.eventlabel}_Meta_Save', save=[])

    print('Finished s03 \n')
    return meta
