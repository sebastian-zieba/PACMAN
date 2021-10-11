import numpy as np
import pysynphot as S 
from ..lib import manageevent as me
#https://pysynphot.readthedocs.io/en/latest/bandpass.html#observation-mode
#https://pysynphot.readthedocs.io/en/latest/appendixb.html#wfc3
#https://pysynphot.readthedocs.io/en/latest/appendixb.html#mjd
import os
from astropy.io import ascii
from scipy.interpolate import interp1d
from ..lib import plots
from tqdm import tqdm


def run03(eventlabel, workdir, meta=None):

    # read in filelist
    filelist_path = meta.workdir + '/filelist.txt'
    if os.path.exists(filelist_path):
        filelist = ascii.read(filelist_path)

    ### Bandpass

    ivisit = filelist['ivisit']
    t_mjd = filelist['t_mjd']

    # stores time when every visit started
    # pysynphot will then download the throughput from the time of observations
    t_mjd_visit_starts = t_mjd[np.unique(ivisit, return_index=True)[1]]

    if meta.grism == 'G141':
        grism = 'g141'
    elif meta.grism == 'G102':
        grism = 'g102'

    if not os.path.exists(meta.workdir + '/ancil/bandpass'):
        os.makedirs(meta.workdir + '/ancil/bandpass')

    for i, t_mjd_visit_start in enumerate(tqdm(t_mjd_visit_starts, desc='Download Bandpass for every visit')):
        bp = S.ObsBandpass('wfc3,ir,{0},mjd#{1}'.format(grism, t_mjd_visit_start))
        wvl, bp_val = bp.binset, bp(bp.binset)
        np.savetxt(meta.workdir + '/ancil/bandpass/bandpass_v{0}.txt'.format(i), list(zip(wvl, bp_val)))

        if meta.save_bandpass_plot or meta.show_bandpass_plot:
            plots.bandpass(wvl, bp_val, grism, i, meta)



    ### Stellar Spectrum

    Teff, logg, MH = meta.Teff, meta.logg, meta.MH

    #def Planck(b_wvl, teff):
    #    h =6.62607004e-34          #m^2/kg/s
    #    c =299792458.0             #m/s
    #    kb =1.38064852e-23         #m^2 kg /s^2 K
    #    return (2.0*h*(c**2)/(b_wvl**5))*(1.0/( np.exp( (h*c)/(b_wvl*kb*teff) )- 1.0))

    if meta.sm == 'phoenix':
        sm = 'phoenix'
    elif meta.sm == 'kurutz':
        sm = 'ck04models'

    # Generate Stellar spectrum
    stellar_spectrum = S.Icat(sm, Teff, MH, logg)
    wvl = stellar_spectrum.wave/1e4 #microns
    flux = stellar_spectrum.flux*1e7/(np.pi)*wvl
    flux = flux/max(flux)

    x = wvl
    y = flux
    f = interp1d(x, y, kind='cubic')


    # Muliply stellar spectrum and bandpass
    for vi in tqdm(np.unique(ivisit), desc='Multiply Bandpass with Stellar Spectrum for every visit'):
        throughput = np.loadtxt(meta.workdir + '/ancil/bandpass/bandpass_v{0}.txt'.format(vi)).T
        wvl_g = throughput[0]/1e4 #microns
        flux_g = throughput[1]
        np.savetxt(meta.workdir + '/ancil/bandpass/refspec_v{0}.txt'.format(vi), list(zip(wvl_g, f(wvl_g)*flux_g)))

        if meta.save_refspec_plot or meta.show_refspec_plot:
            plots.refspec(x, y, wvl_g, flux_g, f, vi, meta)


    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    print('Finished s03 \n')

    return meta
