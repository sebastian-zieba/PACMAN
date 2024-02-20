"""This code reads in the optimally extracted lightcurve and bins it into channels, following Berta '12"""
import time as time_now
from pathlib import Path

import numpy as np
from astropy.io import ascii
from astropy.table import QTable
# from numpy import *
# from pylab import *
from tqdm import tqdm

from .lib import manageevent as me
from .lib import plots
from .lib import sort_nicely as sn
from .lib import util


def run21(eventlabel, workdir: Path, meta=None):
    """This function reads in the lc_spec.txt file with the flux as a
    function of wavelength and bins it into light curves.
    """
    print('Starting s21\n')

    if meta is None:
        meta = me.loadevent(workdir / f'WFC3_{eventlabel}_Meta_Save')

    if meta.use_wvl_list:
        wave_edges = np.array(meta.wvl_edge_list)
        print('Using wvl_edge_list entered by the user: ', wave_edges)
        # now PACMAN can also do overlapping wavelength bins 
        # like in Spake et al., 2018 Nature, Volume 557, Issue 7703, p.68-70 for the bins around 1083 nm
        if len(wave_edges.shape) == 1:
            meta.wvl_bins = int(len(wave_edges)-1)
        elif len(wave_edges.shape) == 2:
            meta.wvl_bins = int(wave_edges.shape[0])
        print('Number of bins:', meta.wvl_bins)
    else:
        meta.wvl_bins = int(meta.wvl_bins)
        wave_edges = np.linspace(meta.wvl_min, meta.wvl_max, meta.wvl_bins+1)*1e4
        print('Number of bins:', meta.wvl_bins)
        print('chosen bin edges:', wave_edges)

    # reads in spectra
    if meta.s21_most_recent_s20:
        # NOTE: If spectra have a specific filetype ending they can also be globbed directly.
        lst_dir = sn.sort_nicely((meta.workdir / "extracted_lc").iterdir())
        # the following line makes sure that only directories starting with a "2" are considered
        # this was implemented after issue #10 was raised (see issue for more info)
        # this works because the dates will always start with a "2"
        lst_dir_new = [lst_dir_i for lst_dir_i in lst_dir if lst_dir_i.name.startswith("2")]
        spec_dir = lst_dir_new[-1]
    else:
        spec_dir = meta.s21_spec_dir_path_s20

    print("Chosen directory with the spectroscopic flux files:", spec_dir)

    # save the mid bin wavelengths into a new file
    table_wvl = QTable(names=('bin', 'wavelengths'))
    if len(wave_edges.shape) == 2:
        wavelengths = np.array([(wave_edges[i][0] + wave_edges[i][1]) / 2. / 1.e4 for i in range(meta.wvl_bins)])
    else:
        wavelengths = np.array([(wave_edges[i] + wave_edges[i+1]) / 2. / 1.e4 for i in range(meta.wvl_bins)])

    d = ascii.read(meta.workdir / "extracted_lc" / spec_dir / "lc_spec.txt")
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

    w_min = w.min()
    w_max = w.max()

    w_hires = np.linspace(w_min, w_max, 10000)
    oversample_factor = len(w_hires)/npix*1.0

    #stores the indices corresponding to the wavelength range in each bin
    wave_inds = []

    if len(wave_edges.shape) == 2:
        wavelengths = np.array([(wave_edges[i][0] + wave_edges[i][1]) / 2. / 1.e4 for i in range(meta.wvl_bins)])
    else:
        wavelengths = np.array([(wave_edges[i] + wave_edges[i+1]) / 2. / 1.e4 for i in range(meta.wvl_bins)])

    if len(wave_edges.shape) == 2:
        for i in range(meta.wvl_bins): 
            wave_inds.append((w_hires >= wave_edges[i][0])&(w_hires <= wave_edges[i][1]))
    else:
        for i in range(meta.wvl_bins): 
            wave_inds.append((w_hires >= wave_edges[i])&(w_hires <= wave_edges[i+1]))

    #for i in range(len(wave_bins)- 1): lo_res_wave_inds.append((w >= wave_bins[i])&(w <= wave_bins[i+1]))

    datetime = time_now.strftime('%Y-%m-%d_%H-%M-%S')
    dirname = meta.workdir / "extracted_sp" / f'bins{meta.wvl_bins}_{datetime}'
    if not dirname.exists():
        dirname.mkdir(parents=True)

    for i in tqdm(range(meta.wvl_bins), desc='***************** Looping over Bins', ascii=True):
        if len(wave_edges.shape) == 2:
            wave = (wave_edges[i][0] + wave_edges[i][1])/2./1.e4
        else:
            wave = (wave_edges[i] + wave_edges[i+1])/2./1.e4
        outname = dirname / f"speclc{wave:.3f}.txt"
        table = QTable(names=('t_mjd', 't_bjd', 't_visit', 't_orbit', 'ivisit', 'iorbit', 'scan', 'spec_opt', 'var_opt', 'wave'))

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

            meanflux, meanerr = util.weighted_mean(fluxes, errs)

            #print(t_mjd, t_bjd, t_visit, t_orbit, ivisit, iorbit, scan, meanflux, meanerr**2, wave, file=outfile)
            #print wave, np.sum(d[j, lo_res_wave_inds[i],2])
            table.add_row([t_mjd_i, t_bjd_i, t_visit_i, t_orbit_i, ivisit_i, iorbit_i, scan_i, meanflux, meanerr**2, wave])

    #print wave, 1.0*sum(wave_inds)/len(w_hires), meanflux, meanerr
        ascii.write(table, outname, format='ecsv', overwrite=True)

    print(f'Saved light curve(s) in {dirname}')
    plots.plot_wvl_bins(w_hires, f_interp, wave_edges, meta.wvl_bins, dirname)

    print('Saving Wavelength bin file')
    for idx, wavelengths_i in enumerate(wavelengths):
        table_wvl.add_row([idx, wavelengths_i])
    ascii.write(table_wvl, dirname / 'wvl_table.dat', format='rst', overwrite=True)

    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir / f'WFC3_{meta.eventlabel}_Meta_Save', save=[])

    print('Finished s21 \n')
    return meta
