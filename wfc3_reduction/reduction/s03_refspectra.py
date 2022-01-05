import matplotlib.pyplot as plt
import numpy as np
import pysynphot as S
from ..lib import manageevent as me
# https://pysynphot.readthedocs.io/en/latest/bandpass.html#observation-mode
# https://pysynphot.readthedocs.io/en/latest/appendixb.html#wfc3
# https://pysynphot.readthedocs.io/en/latest/appendixb.html#mjd
import os
from astropy.io import ascii
from scipy.interpolate import interp1d
from ..lib import plots
from tqdm import tqdm
from ..lib import stellar_spectrum


def binning(x_input, y_input, x_ref):
    """
    A binning function which bins x_input and y_input values and uses x_ref to determine boarders.

    TODO improve the function
    TODO check if oversampling is necessary
    """
    #print(x_input)
    n = 4
    x_input_mids = []
    y_input_mids = []
    for i in range(len(x_input)-1):
        x = x_input[i]
        y = x_input[i+1]
        step = (y - x) / (n - 1)
        x_input_mid = [x + step * i for i in np.arange(1,n-1)]
        x_input_mids.append(x_input_mid)

        x = y_input[i]
        y = y_input[i+1]
        step = (y - x) / (n - 1)
        y_input_mid = [x + step * i for i in np.arange(1,n-1)]
        y_input_mids.append(y_input_mid)

    def flatten(t):
        return np.array([item for sublist in t for item in sublist])

    x_input_mids = flatten(x_input_mids)
    y_input_mids = flatten(y_input_mids)

    x_input = np.concatenate((x_input, x_input_mids))
    y_input = np.concatenate((y_input, y_input_mids))
    x_input_sort = np.argsort(x_input)
    x_input = x_input[x_input_sort]
    y_input = y_input[x_input_sort]

    x_ref0 = np.array([x_ref[0] - (x_ref[1] - x_ref[0])])
    x_ref1 = np.array([x_ref[-1] - (x_ref[-2] - x_ref[-1])])
    x_ref_new = np.concatenate((x_ref0, x_ref, x_ref1))
    x_ref_edges = (x_ref_new[:-1] + x_ref_new[1:]) / 2
    x_input_cuts = []
    y_input_cuts = []
    for i in range(len(x_ref_edges) - 1):
        x_input_cut = []
        y_input_cut = []
        for j in range(len(x_input)):
            if x_ref_edges[i] <= x_input[j] < x_ref_edges[i + 1]:
                #print(x_ref_edges[i], x_ref_edges[i+1])
                x_input_cut.append(x_input[j])
                y_input_cut.append(y_input[j])
        x_input_cuts.append(np.array(x_input_cut))
        y_input_cuts.append(np.array(y_input_cut))
    x_input_cuts = np.array(x_input_cuts, dtype=object)
    y_input_cuts = np.array(y_input_cuts, dtype=object)
    #print([len(ii) for ii in x_input_cuts])
    #for ii in range(len(x_input_cuts)):
        #if len(x_input_cuts[ii]) == 0:
            #print(ii, x_input_cuts[ii])
    x_input_means = np.array([np.mean(ii) for ii in x_input_cuts])
    y_input_means = np.array([np.mean(ii) for ii in y_input_cuts])
    return x_input_means, y_input_means


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
	   meta object with all the meta data stored in s01

	History
	-------
	Written by Sebastian Zieba      December 2021
    """

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    if meta.grism == 'G141':
        grism = 'g141'
    elif meta.grism == 'G102':
        grism = 'g102'

    bp_wvl, bp_val = np.loadtxt(meta.ancildir + '/bandpass/bandpass_{0}.txt'.format(grism)).T

    # if meta.save_bandpass_plot or meta.show_bandpass_plot:
    #    plots.bandpass(bp_wvl, bp_val, grism, 0, meta)

    ### Stellar Spectrum
    Teff, logg, MH = meta.Teff, meta.logg, meta.MH
    meta.sm='ck04models'
    print(meta.sm)
    print(meta.sm in ['k93models', 'ck04models', 'phoenix'])
    if meta.sm in ['k93models', 'ck04models', 'phoenix']:
        sm_wvl, sm_flux = stellar_spectrum.get_sm(meta, MH, logg, Teff)
    elif meta.sm == 'blackbody':
        sm_wvl, sm_flux = stellar_spectrum.get_bb(Teff)
    else:
        print('You have not entered a valid stellar model.\n'
              'Options are k93models, ck04models, phoenix or blackbody \n'
              'We proceed with the bb')
        sm_wvl, sm_flux = stellar_spectrum.get_bb(Teff)

    bp_wvl = bp_wvl * 1e-10
    bp_val = bp_val / max(bp_val)
    sm_flux = sm_flux #* sm_wvl  # in order to convert from W/m^3/sr units to W/m^2/sr
    sm_flux = sm_flux / max(sm_flux)

    plt.plot(bp_wvl*1e6, bp_val)
    plt.xlim(0.6, 1.7)
    plt.savefig(meta.workdir + '/ancil/bandpass/bandpass.png')
    plt.close()
    plt.plot(sm_wvl*1e6, sm_flux)
    plt.xlim(0.6, 1.7)
    plt.savefig(meta.workdir + '/ancil/bandpass/sm.png')
    plt.close()

    sm_wvl_binned, sm_flux_binned = binning(sm_wvl, sm_flux, bp_wvl)

    # sm_wvl_cut = []
    # sm_flux_cut = []
    # for i, sm_wvl_i in enumerate(sm_wvl_binned):
    #     if bp_wvl[0] <= sm_wvl_i <= bp_wvl[-1]:
    #         sm_wvl_cut.append(sm_wvl_binned[i])
    #         sm_flux_cut.append(sm_flux_binned[i])

    sm_wvl_cut = np.array(sm_wvl_binned)
    sm_flux_cut = np.array(sm_flux_binned)

    f = interp1d(sm_wvl_cut, sm_flux_cut, kind='linear')

    #print(sm_wvl_cut[0], sm_wvl_cut[-1])
    #print(bp_wvl[0], bp_wvl[-1])

    wvl_ref = bp_wvl[1:-1]
    flux_ref = f(bp_wvl[1:-1]) * bp_val[1:-1]
    flux_ref = flux_ref / max(flux_ref)

    if not os.path.exists(meta.workdir + '/ancil/bandpass/'):
        os.mkdir(meta.workdir + '/ancil/bandpass/')
    np.savetxt(meta.workdir + '/ancil/bandpass/refspec.txt', list(zip(wvl_ref, flux_ref)))

    if meta.save_refspec_plot or meta.show_refspec_plot:
        plots.refspec(bp_wvl, bp_val, sm_wvl, sm_flux, wvl_ref, flux_ref, meta)

    #
    # # Generate Stellar spectrum
    # stellar_spectrum = S.Icat(sm, Teff, MH, logg)
    # wvl = stellar_spectrum.wave/1e4 #microns
    # flux = stellar_spectrum.flux*1e7/(np.pi)*wvl
    # flux = flux/max(flux)
    #
    # x = wvl
    # y = flux
    # f = interp1d(x, y, kind='cubic')
    #
    #
    # # Muliply stellar spectrum and bandpass
    # for vi in tqdm(np.unique(ivisit), desc='Multiply Bandpass with Stellar Spectrum for every visit'):
    #     throughput = np.loadtxt(meta.workdir + '/ancil/bandpass/bandpass_v{0}.txt'.format(vi)).T
    #     wvl_g = throughput[0]/1e4 #microns
    #     flux_g = throughput[1]
    #     np.savetxt(meta.workdir + '/ancil/bandpass/refspec_v{0}.txt'.format(vi), list(zip(wvl_g, f(wvl_g)*flux_g)))
    #
    #     if meta.save_refspec_plot or meta.show_refspec_plot:
    #         plots.refspec(x, y, wvl_g, flux_g, f, vi, meta)
    #
    #
    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    print('Finished s03 \n')

    return meta
