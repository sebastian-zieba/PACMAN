#This code reads in the optimally extracted lightcurve and bins it into channels 5 pixels wide, following Berta '12
import numpy as np
#from numpy import *
#from pylab import *
from astropy.io import ascii
from scipy import signal
import os
import time as time_now
from astropy.table import QTable
from tqdm import tqdm
from .lib import plots
from .lib import sort_nicely as sn
from .lib import manageevent as me
from astropy.table import QTable


def run21(eventlabel, workdir, meta=None):
    """
    This function reads in the lc_spec.txt file with the flux as a funtion of wavelength and bins it into light curves.
    """
    print('Starting s21\n')

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    def weighted_mean(data, err):				#calculates the weighted mean for data points data with std err
        weights = 1.0/err**2.
        mu = np.sum(data*weights)/np.sum(weights)
        var = 1.0/np.sum(weights)
        return [mu, np.sqrt(var)]				#returns weighted mean and variance


    if meta.use_wvl_list:
        print(meta.wvl_edge_list)
        wave_edges = np.array(meta.wvl_edge_list)
        meta.wvl_bins = int(len(wave_edges)-1)
        print('Number of bins:', meta.wvl_bins)
    else:
        meta.wvl_bins = int(meta.wvl_bins)
        wave_edges = np.linspace(meta.wvl_min, meta.wvl_max, meta.wvl_bins+1)*1e4
        print('Number of bins:', meta.wvl_bins)
        print('chosen bin edges:', wave_edges)

    #reads in spectra
    if meta.s21_most_recent_s20:
        lst_dir = os.listdir(meta.workdir + "/extracted_lc/")
        lst_dir = sn.sort_nicely(lst_dir)
        spec_dir = lst_dir[-1]
    else:
        spec_dir = meta.s21_spec_dir_path_s20

    print("Chosen directory with the spectroscopic flux files:", spec_dir)

    # save the mid bin wavelengths into a new file
    table_wvl = QTable(names=('bin', 'wavelengths'))
    wavelengths = np.array([(wave_edges[i] + wave_edges[i+1]) / 2. / 1.e4 for i in range(len(wave_edges) - 1)])

    d = ascii.read(meta.workdir + "/extracted_lc/" + spec_dir + "/lc_spec.txt")
    d = np.array([d[i].data for i in d.colnames])

    nexp = meta.nexp		            #number of exposures
    npix = meta.npix#meta.BEAMA_f - meta.BEAMA_i  #width of spectrum in pixels (BEAMA_f - BEAMA_i)
    #d = d.reshape(nexp , npix,  -1)			#reshapes array by exposure

    t_mjd, t_bjd = d[0].reshape(nexp, npix), d[1].reshape(nexp, npix)
    t_visit, t_orbit = d[2].reshape(nexp, npix), d[3].reshape(nexp, npix)
    ivisit, iorbit = d[4].reshape(nexp, npix), d[5].reshape(nexp, npix)
    scan = d[6].reshape(nexp, npix)
    spec_opt, var_opt = d[7].reshape(nexp, npix), d[8].reshape(nexp, npix)
    w = d[9].reshape(nexp, npix) # d[0,:, 4]
    #print(w[0])
    #f = d[0, :, 2]

    w_min = w.min()#max(w[:,0])
    w_max = w.max()#min(w[:,-1])
    #print(w_min, w_max)
    #print(w.min(), w.max())
    #w_hires = np.linspace(w.min(), w.max(), 10000)
    w_hires = np.linspace(w_min, w_max, 10000)
    oversample_factor = len(w_hires)/npix*1.0
    #print(oversample_factor)
    #stores the indices corresponding to the wavelength range in each bin
    wave_inds = []
    #lo_res_wave_inds = []
    for i in range(len(wave_edges)- 1): wave_inds.append((w_hires >= wave_edges[i])&(w_hires <= wave_edges[i+1]))
    #for i in range(len(wave_bins)- 1): lo_res_wave_inds.append((w >= wave_bins[i])&(w <= wave_bins[i+1]))

    datetime = time_now.strftime('%Y-%m-%d_%H-%M-%S')
    dirname = meta.workdir + "/extracted_sp/" + 'bins{0}_'.format(meta.wvl_bins) + datetime
    if not os.path.exists(dirname): os.makedirs(dirname)

    for i in tqdm(range(len(wave_edges) - 1), desc='***************** Looping over Bins', ascii=True):

        wave = (wave_edges[i] + wave_edges[i+1])/2./1.e4
        outname = dirname + "/speclc" + "{0:.3f}".format(wave)+".txt"
        #outname = "wasp33b_" + "{0:.4f}".format(wave)+".txt"
        #outfile = open(outname, 'w')
        #print(sum(wave_inds[i]))
        table = QTable(names=('t_mjd', 't_bjd', 't_visit', 't_orbit', 'ivisit', 'iorbit', 'scan', 'spec_opt', 'var_opt', 'wave'))

        #print('#t_mjd', '\t', 't_bjd', '\t', 't_visit', '\t', 't_orbit', '\t', 'ivisit', '\t', 'iorbit', '\t', 'scan', '\t', 'spec_opt', '\t', 'var_opt', '\t','wave', file=outfile)
        for j in range(nexp):
            t_mjd_i, t_bjd_i = t_mjd[j][0], t_bjd[j][0]
            t_visit_i, t_orbit_i = t_visit[j][0], t_orbit[j][0]
            ivisit_i, iorbit_i = ivisit[j][0], iorbit[j][0]
            scan_i = scan[j][0]
            spec_opt_i,  var_opt_i = spec_opt[j], var_opt[j]
            w_i = w[j]

            f_interp = np.interp(w_hires, w_i, spec_opt_i)
            variance_interp = np.interp(w_hires, w_i, var_opt_i)

            #accounts for decrease in precision when spectrum is oversampled
            variance_interp *= oversample_factor

            fluxes = f_interp[wave_inds[i]]
            errs = np.sqrt(variance_interp[wave_inds[i]])

            meanflux, meanerr = weighted_mean(fluxes, errs)

            #print(t_mjd, t_bjd, t_visit, t_orbit, ivisit, iorbit, scan, meanflux, meanerr**2, wave, file=outfile)
            #print wave, np.sum(d[j, lo_res_wave_inds[i],2])
            table.add_row([t_mjd_i, t_bjd_i, t_visit_i, t_orbit_i, ivisit_i, iorbit_i, scan_i, meanflux, meanerr**2, wave])

    #print wave, 1.0*sum(wave_inds)/len(w_hires), meanflux, meanerr
        ascii.write(table, outname, format='ecsv', overwrite=True)

    print('Saved light curve(s) in {0}'.format(dirname))

    plots.plot_wvl_bins(w_hires, f_interp, wave_edges, meta.wvl_bins, dirname)

    print('Saving Wavelength bin file')
    for idx, wavelengths_i in enumerate(wavelengths):
        table_wvl.add_row([idx, wavelengths_i])
    ascii.write(table_wvl, dirname + '/wvl_table.dat', format='rst', overwrite=True)

    print('Saving Metadata')
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    print('Finished s21 \n')

    return meta
