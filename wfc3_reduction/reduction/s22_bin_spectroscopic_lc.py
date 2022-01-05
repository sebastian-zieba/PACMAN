#This code reads in the optimally extracted lightcurve and bins it into channels 5 pixels wide, following Berta '12
import numpy as np
from numpy import *
from pylab import *
from astropy.io import ascii
from scipy import signal
import matplotlib.pyplot as plt
import os
from ..lib import manageevent as me
from ..lib import util
import time as time_now
from astropy.table import QTable
from ..lib import plots


def run22(eventlabel, workdir, meta=None):

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    def weighted_mean(data, err):				#calculates the weighted mean for data points data with std err
        weights = 1.0/err**2.
        mu = np.sum(data*weights)/np.sum(weights)
        var = 1.0/np.sum(weights)
        return [mu, np.sqrt(var)]				#returns weighted mean and variance


    for wvl_bins in meta.wvl_bins:
        print(meta.wvl_bins)
        #what bins do you want?
        #wave_bins = np.linspace(1.125, 1.65, 22)*1e4
        wave_bins = np.linspace(meta.wvl_min, meta.wvl_max, wvl_bins)*1e4
        print(wave_bins)

        #file_mario = '/home/zieba/Documents/L98-59_bins_2.txt'
        #wvl_mario0, wvl_mario1,_,_,_,_=np.loadtxt(file_mario).T
        #wave_bins = np.concatenate((wvl_mario0, [wvl_mario1[-1]]))+154
        #a = 8000
        #b = 10809
        #c = 10858
        #d = 11400
        #ab = np.linspace(a, b, 13)
        #cd = np.linspace(c, d, 4)

        #wave_bins = np.concatenate((ab,cd))
        print(wave_bins)

        #wave_bins = np.linspace(1.139, 1.631, 12)*1e4
        #wave_bins = np.array([11400, 11800, 12200,12600,13000, 13400, 13800, 14200, 14600, 15000, 15400, 15800, 16200])

        #Tom Evan wave bins
        #wave_bins = np.array([1.129595435219961752e+00, 1.180470435219961756e+00, 1.231345435219961759e+00, 1.282220435219961763e+00, 1.333095435219961766e+00, 1.383970435219961770e+00, 1.434845435219961773e+00, 1.485720435219961777e+00, 1.536595435219961781e+00, 1.587470435219961784e+00, 1.638345435219961788e+00])*1.e4

        #reads in spectra

        lst = os.listdir(meta.workdir + "/extracted_lc/")
        print(lst)
        print(lst.sort())
        d = ascii.read(meta.workdir + "/extracted_lc/" + lst[0] + "/lc_spec.txt")
        d = np.array([d[i].data for i in d.colnames])

        nexp = meta.nexp		#number of exposures
        npix = meta.BEAMA_f - meta.BEAMA_i #181  #width of spectrum in pixels (BEAMA_f - BEAMA_i)
        #d = d.reshape(nexp , npix,  -1)			#reshapes array by exposure

        t_mjd, t_bjd = d[0].reshape(nexp, npix), d[1].reshape(nexp, npix)
        t_visit, t_orbit = d[2].reshape(nexp, npix), d[3].reshape(nexp, npix)
        ivisit, iorbit = d[4].reshape(nexp, npix), d[5].reshape(nexp, npix)
        scan = d[6].reshape(nexp, npix)
        spec_opt, var_opt = d[7].reshape(nexp, npix), d[8].reshape(nexp, npix)
        w = d[9].reshape(nexp, npix) # d[0,:, 4]
        #f = d[0, :, 2]

        w_min = max(w[:,0])
        w_max = min(w[:,-1])

        #w_hires = np.linspace(w.min(), w.max(), 10000)
        w_hires = np.linspace(w_min, w_max, 10000)
        oversample_factor = len(w_hires)/npix*1.0
        print(oversample_factor)
        #stores the indices corresponding to the wavelength range in each bin
        wave_inds = []
        #lo_res_wave_inds = []
        for i in range(len(wave_bins)- 1): wave_inds.append((w_hires >= wave_bins[i])&(w_hires <= wave_bins[i+1]))
        #for i in range(len(wave_bins)- 1): lo_res_wave_inds.append((w >= wave_bins[i])&(w <= wave_bins[i+1]))

        datetime = time_now.strftime('%Y-%m-%d_%H-%M-%S')
        dirname = meta.workdir + "/extracted_sp/" + 'bins{0}_'.format(wvl_bins) + datetime
        if not os.path.exists(dirname): os.makedirs(dirname)

        for i in range(len(wave_bins) - 1):
            wave = (wave_bins[i] + wave_bins[i+1])/2./1.e4
            outname = dirname + "/speclc" + "{0:.3f}".format(wave)+".txt"
            #outname = "wasp33b_" + "{0:.4f}".format(wave)+".txt"
            #outfile = open(outname, 'w')

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

        plots.plot_wvl_bins(w_hires, f_interp, wave_bins, wvl_bins, dirname)

        #outfile.close()

    print('Finished s22 \n')

    return meta