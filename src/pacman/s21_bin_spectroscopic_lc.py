import time
import shutil
from pathlib import Path

import numpy as np
from astropy.io import ascii
from astropy.table import QTable
from tqdm import tqdm

from .lib import manageevent as me
from .lib import plots
from .lib import util
from .lib import logedit
from .lib import read_pcf as rd


def run21(pcf_path: Path, meta=None):
    """
    This function reads in the lc_spec.txt file with the flux as a
    function of wavelength and bins it into light curves.
    """

    meta, log = util.setup_stage(
        pcf_path=pcf_path,
        stage_num="21",
        previous_stage_num="20",
        copy_filelist=True,
        copy_xrefyref=True,
        copy_ancil=True,
        copy_extracted_lc=True,
        meta=meta,
    )

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
    spec_dir = meta.workdir / "extracted_lc"
    log.writelog(f"Using spectroscopic flux files from: {spec_dir}")

    print("Chosen directory with the spectroscopic flux files:", spec_dir)

    d = ascii.read(str(spec_dir / "lc_spec.txt"))
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

    datetime = time.strftime("%Y-%m-%d_%H-%M-%S")
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

        plots.light_curve_errorbar(
            outname,
            meta.workdir / "figs" / "s21_lightcurves",
            f"speclc{wave:.3f}.png",
            title=f"Spectroscopic light curve: {wave:.3f} micron",
        )

    log.writelog(f"Saved light curve(s) in {dirname}")

    fig_dir = meta.workdir / "figs" / "s21_lightcurves"
    fig_dir.mkdir(parents=True, exist_ok=True)
    plots.plot_wvl_bins(w_hires, f_interp, wave_edges, meta.wvl_bins, fig_dir)

    # save the mid bin wavelengths into a new file
    log.writelog("Saving Wavelength bin file")
    table_wvl = QTable(
        names=("bin", "wavelength", "half_width", "lower_edge", "upper_edge")
    )
    for idx in range(meta.wvl_bins):
        if len(wave_edges.shape) == 2:
            lower_edge = wave_edges[idx][0] / 1.0e4
            upper_edge = wave_edges[idx][1] / 1.0e4
        else:
            lower_edge = wave_edges[idx] / 1.0e4
            upper_edge = wave_edges[idx + 1] / 1.0e4

        wavelength = 0.5 * (lower_edge + upper_edge)
        half_width = 0.5 * (upper_edge - lower_edge)

        table_wvl.add_row([idx, wavelength, half_width, lower_edge, upper_edge])

    ascii.write(table_wvl, dirname / "wvl_table.dat", format="rst", overwrite=True)

    # Save results
    log.writelog("Saving Metadata")
    me.saveevent(meta, meta.workdir / "WFC3_Meta_Save", save=[])

    log.writelog("Finished s21 \n")
    log.closelog()
    return meta
