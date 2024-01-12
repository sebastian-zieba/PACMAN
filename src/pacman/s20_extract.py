import sys
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
from astropy.io import ascii, fits
from astropy.table import QTable
from tqdm import tqdm

from .lib import manageevent as me
from .lib import optextr
from .lib import plots
from .lib import util


def run20(eventlabel, workdir: Path, meta=None):
    """This function extracts the spectrum and saves the total flux and the
    flux as a function of wavelength into files.
    """
    print('Starting s20')

    if meta is None:
        meta = me.loadevent(workdir / f'WFC3_{eventlabel}_Meta_Save')

    # NOTE: Load in more information into meta
    meta = util.ancil(meta, s20=True)

    ###############################################################################################################################################################
    # NOTE: STEP 0: Set up files and directories
    ###############################################################################################################################################################

    dirname = meta.workdir / "extracted_lc" /\
        datetime.strftime(datetime.now(), '%Y-%m-%d_%H-%M-%S')
    if not dirname.exists():
        dirname.mkdir(parents=True)

    # NOTE: Let's have a copy of the pcf in the extracted_lc directory
    # This copy is just for the user to know what parameters they used when
    # running s20
    shutil.copy(meta.workdir / 'obs_par.pcf', dirname)

    # NOTE: Initialize the astropy tables where we will save the extracted
    # spectra
    if meta.output:
        table_white = QTable(names=('t_mjd', 't_bjd', 't_visit','t_orbit', 'ivisit', 'iorbit', 'scan', 'spec_opt', 'var_opt', 'spec_box', 'var_box'))
        table_spec = QTable(names=('t_mjd', 't_bjd', 't_visit','t_orbit', 'ivisit', 'iorbit', 'scan', 'spec_opt', 'var_opt', 'template_waves'))
        table_diagnostics = QTable(names=('nspectra', 't_mjd', 'numoutliers', 'skymedian', "# nans"))

    files_sp = meta.files_sp     # spectra files
    nspectra = 0                 # iterator variable to track number of spectra reduced

    # NOTE: Only do the first N files, if wanted by the user
    if meta.s20_testing:
        meta.nexp = meta.n_testing
        print('Running s20 in testing mode...')
    else:
        meta.nexp = len(files_sp)    # number of exposures

    # TODO: add pcf flag to only save the plots for the very first file

    # NOTE: The following lists are used for diagnostic plots
    if meta.save_utr_aper_evo_plot or meta.show_utr_aper_evo_plot:
        peaks_all = []
    if meta.save_bkg_evo_plot or meta.show_bkg_evo_plot:
        bkg_evo = []
    if meta.save_sp1d_diff_plot or meta.show_sp1d_diff_plot:
        sp1d_all = []
        wvl_hires = np.linspace(7000, 17800, 1000)
    if meta.save_drift_plot or meta.show_drift_plot:
        leastsq_res_all = []

    print('in total #visits, #orbits:', (meta.nvisit, meta.norbit), '\n')
    if meta.background_box:
        print('Using Background Box')

    # NOTE: In order to have the correct order of print() with tqdm,
    # i added file=sys.stdout
    # source: https://stackoverflow.com/questions/36986929/redirect-print-command-in-python-script-through-tqdm-write
    for i in tqdm(np.arange(meta.nexp, dtype=int), desc='***************** Looping over files', file=sys.stdout):#tqdm(np.arange(len(files_sp), dtype=int)):
        f = files_sp[i]             # Current file
        print(f"\nFilename: {f}")
        d = fits.open(f)            # Opens the file
        scan = meta.scans_sp[i]     # Scan direction of the spectrum

        # NOTE: Plot with good visible background
        if meta.save_sp2d_plot or meta.show_sp2d_plot:
            plots.sp2d(d, meta, i)

        visnum, orbnum = meta.ivisit_sp[i], meta.iorbit_sp_cumulative[i]     # Current visit and cumulative orbit number
        print('current visit, orbit: ', (visnum, orbnum))

        # NOTE: Plot trace y pos in the trace plot is the position of the DI
        if meta.save_trace_plot or meta.show_trace_plot:
            plots.trace(d, meta, visnum, orbnum, i)

        # TODO: SPEED UP: calculation of the start and end of the trace could be moved to util.py. It's also used in plots.plot_trace. Also in plots.utr
        cmin = int(meta.refpix[orbnum, 2] + meta.POSTARG1/meta.platescale) + meta.BEAMA_i + meta.LTV1                      #determines left column for extraction (beginning of the trace)
        cmax = min(int(meta.refpix[orbnum, 2] + meta.POSTARG1/meta.platescale) + meta.BEAMA_f + meta.LTV1, meta.subarray_size-5)     #right column (end of trace, or edge of detector)
        rmin, rmax = int(meta.rmin), int(meta.rmax)     # Top and bottom row for extraction (specified in obs_par.txt)
        meta.npix = cmax - cmin

        # TODO: SPEED UP: This calculation is done for every file again from new!
        # TODO: Move it to util.py so its done just once
        meta.wave_grid = util.get_wave_grid(meta)   # Gets grid of wavelength solutions for each orbit and row

        M = np.ones_like(d[1].data[rmin:rmax, cmin:cmax])                        #mask for bad pixels
        flatfield = util.get_flatfield(meta)
        bpix = d[3].data[rmin:rmax,cmin:cmax]
        badpixind =  (bpix==4)|(bpix==512)|(flatfield[orbnum][rmin:rmax, cmin:cmax] == -1.)    #selects bad pixels
        # print('bad pixels', sum(bpix==4), sum(bpix==512),sum(flatfield[orbnum][rmin:rmax, cmin:cmax] == -1.), sum(badpixind))
        if i == 0: plots.badmask_2d((bpix==4), (bpix==512), (flatfield[orbnum][rmin:rmax, cmin:cmax] == -1), meta, i)
        M[badpixind] = 0.0                                        #initializes bad pixel mask
        # NOTE: Store number of bad pixels
        spec_box = np.zeros(cmax - cmin)                                #box extracted standard spectrum
        spec_opt = np.zeros(cmax - cmin)                                #optimally extracted spectrum
        var_box = np.zeros(cmax - cmin)                              #box spectrum variance
        var_opt = np.zeros(cmax - cmin)                                #optimal spectrum variance

        #########################################################################################################################################################
        # NOTE: Loops over up-the-ramp-samples (skipping first two very short exposures); gets all needed input for optextr routine                    #
        #########################################################################################################################################################
        # In order to not print a new line with tqdm every time, I added
        # leave=True, position=0
        # source: https://stackoverflow.com/questions/41707229/tqdm-printing-to-newline
        if len(set(meta.scans_sp)) == 1:
            for ii in tqdm(np.arange(d[0].header['nsamp']-1, dtype=int), desc='--- Looping over up-the-ramp-samples',
                           leave=True, position=0):

                fullframe = d[ii * 5 + 1].data
                frame = fullframe[rmin:rmax, cmin:cmax]

                # NOTE: Calculate aperture
                peaks = util.peak_finder(frame, i, ii, orbnum, meta)

                # NOTE: Stores the locations of the peaks for every file and
                # up-the-ramp-samples
                if meta.save_utr_aper_evo_plot or meta.show_utr_aper_evo_plot:
                    peaks_all.append(peaks)

                ### BACKGROUND SUBTRACTION
                if not meta.background_box:
                    below_threshold = fullframe < meta.background_thld # mask with all pixels below the user defined threshold
                    skymedian = np.median(fullframe[below_threshold].flatten())  # estimates the background counts by taking the flux median of the pixels below the flux threshold
                    if meta.save_bkg_evo_plot or meta.show_bkg_evo_plot:
                        bkg_evo.append(skymedian)
                    skyvar = util.median_abs_dev(fullframe[below_threshold].flatten())  # variance for the background count estimate
                    if meta.save_bkg_hist_plot or meta.show_bkg_hist_plot:
                        plots.bkg_hist(fullframe, skymedian, meta, i, ii)
                else:
                    bg_rmin, bg_rmax, bg_cmin, bg_cmax = int(meta.bg_rmin), int(meta.bg_rmax), int(meta.bg_cmin), int(meta.bg_cmax),
                    skymedian = np.median(fullframe[bg_rmin:bg_rmax, bg_cmin:bg_cmax].flatten())
                    if meta.save_bkg_evo_plot or meta.show_bkg_evo_plot:
                        bkg_evo.append(skymedian)
                    skyvar = util.median_abs_dev(fullframe[bg_rmin:bg_rmax, bg_cmin:bg_cmax].flatten())

                frame = frame - skymedian 

                spectrum = frame[max(min(peaks) - meta.window, 0):min(max(peaks) + meta.window, rmax), :]

                err = np.zeros_like(spectrum) + float(meta.rdnoise) ** 2 + skyvar
                var = abs(spectrum) + float(
                    meta.rdnoise) ** 2 + skyvar     # Variance estimate: Poisson noise from photon counts (first term)  + readnoise (factor of 2 for differencing) + skyvar
                spec_box_0 = spectrum.sum(axis=0)   # Initial box-extracted spectrum
                var_box_0 = var.sum(axis=0)         # Initial variance guess
                # Mnew = np.ones_like(M[max(min(peaks) - meta.window, 0):min(max(peaks) + meta.window, rmax),:])
                # Mnew = M[max(peaks_mid - 110, 0):min(peaks_mid + 110, rmax),:]
                Mnew = M[max(min(peaks) - meta.window, 0):min(max(peaks) + meta.window, rmax), :]
                # TODO: Just use meta to reduce the number of parameters which are given to optextr
                if meta.opt_extract == True:
                    [f_opt_0, var_opt_0, numoutliers] = optextr.optextr(spectrum, err, spec_box_0, var_box_0, Mnew,
                                                                        meta.nsmooth, meta.sig_cut,
                                                                        meta.save_optextr_plot, i, ii, meta)
                else:
                    [f_opt, var_opt] = [spec_box_0, var_box_0]

                # NOTE: Sums up spectra and variance for all the differenced
                # images
                spec_opt += f_opt_0
                var_opt += var_opt_0
                spec_box += spec_box_0
                var_box += var_box_0

        elif len(set(meta.scans_sp)) == 2:
            # NOTE: If needed: Calculate spectrum shift in spatial direction (rows) with respect to the first exposure in the data-set (i.e. first visit, first orbit, first exposure)
            if meta.calculate_rowshift:
                #get spatial profile of exposure of first readout
                row_spat, flux_spat = util.get_boxspat(meta,i,d,rmin,rmax,cmin,cmax,(d[0].header['nsamp']-1)*5+1) #row_spat is always np.linspace(rmin,rmax+1,(rmax-rmin)*1000)
                # NOTE: determine rowshift and store it in meta
                if i == 0:
                    #initialize arrays
                    meta.refprofile_rowshift = np.zeros((2,len(row_spat)))
                    meta.rowshift = np.zeros(meta.nexp)

                if i < 2:
                    #store reference profiles of the first exposures
                    meta.refprofile_rowshift[meta.scans_sp[i]] = flux_spat    
                else:
                    #fit spatial profiles to reference profiles
                    rowshift_exp = util.calculate_rowshift(meta,row_spat, meta.refprofile_rowshift[meta.scans_sp[i]],row_spat[10000:-10000], flux_spat[10000:-10000],i)
                    meta.rowshift[i] = rowshift_exp
                
                # NOTE: calculate stretch of last readout with respect to first exposure in the data-set 
                if meta.save_rowshift_stretch_plot or meta.show_rowshift_stretch_plot:
                    row_spat, flux_spat = util.get_boxspat(meta,i,d,rmin,rmax,cmin,cmax,1) #row_spat is always np.linspace(rmin,rmax+1,(rmax-rmin)*1000)
                    # NOTE: determine rowshift and store it in meta
                    if i == 0:
                        #initialize arrays
                        meta.refprofile_stretch = np.zeros((2,len(row_spat)))
                        meta.stretch = np.ones(meta.nexp)

                    if i < 2:
                        #store reference profiles of the first exposures
                        meta.refprofile_stretch[meta.scans_sp[i]] = flux_spat    
                    else:
                        #fit spatial profiles to reference profiles
                        stretch_exp = util.calculate_stretch(meta,row_spat, meta.refprofile_stretch[meta.scans_sp[i]],row_spat[10000:-10000], flux_spat[10000:-10000],i)  
                        meta.stretch[i] = stretch_exp

            for ii in tqdm(np.arange(d[0].header['nsamp']-2, dtype=int), desc='--- Looping over up-the-ramp-samples', leave=True, position=0):
                diff = d[ii*5 + 1].data[rmin:rmax,cmin:cmax] - d[ii*5 + 6].data[rmin:rmax,cmin:cmax]    # Creates image that is the difference between successive scans

                # NOTE: Calculate aperture
                peaks = util.peak_finder(diff, i, ii, orbnum, meta)

                # NOTE: Stores the locations of the peaks for every file
                # and up-the-ramp-samples
                if meta.save_utr_aper_evo_plot or meta.show_utr_aper_evo_plot:
                    peaks_all.append(peaks)

                # NOTE: Sstimates sky background and variance
                fullframe_diff = d[ii*5 + 1].data - d[ii*5 + 6].data    # Fullframe difference between successive scans

                # NOTE: Background subtraction
                if not meta.background_box:
                    below_threshold = fullframe_diff < meta.background_thld             # Mask with all pixels below the user defined threshold
                    skymedian = np.median(fullframe_diff[below_threshold].flatten())    # Estimates the background counts by taking the flux median of the pixels below the flux threshold
                    if meta.save_bkg_evo_plot or meta.show_bkg_evo_plot:
                        bkg_evo.append(skymedian)
                    skyvar = util.median_abs_dev(fullframe_diff[below_threshold].flatten())  # Variance for the background count estimate
                    if meta.save_bkg_hist_plot or meta.show_bkg_hist_plot:
                        plots.bkg_hist(fullframe_diff, skymedian, meta, i, ii)
                else:
                    bg_rmin, bg_rmax, bg_cmin, bg_cmax = int(meta.bg_rmin), int(meta.bg_rmax), int(meta.bg_cmin), int(meta.bg_cmax),
                    skymedian = np.median(fullframe_diff[bg_rmin:bg_rmax, bg_cmin:bg_cmax].flatten())
                    if meta.save_bkg_evo_plot or meta.show_bkg_evo_plot:
                        bkg_evo.append(skymedian)
                    skyvar = util.median_abs_dev(fullframe_diff[bg_rmin:bg_rmax, bg_cmin:bg_cmax].flatten())

                diff = diff - skymedian                                    #subtracts the background

                # print(peaks[0]-peaks[1])

                # peaks_mid = int((peaks[0]+peaks[1])/2)
                # NOTE: Selects postage stamp centered around spectrum
                # we use a bit more data by using the user defined window
                # spectrum = diff[max(peaks_mid - 110, 0):min(peaks_mid + 110, rmax),:]
                spectrum = diff[max(min(peaks) - meta.window, 0):min(max(peaks) + meta.window, rmax),:]

                err = np.zeros_like(spectrum) + float(meta.rdnoise)**2 + skyvar
                var = abs(spectrum) + float(meta.rdnoise)**2 + skyvar          # variance estimate: Poisson noise from photon counts (first term)  + readnoise (factor of 2 for differencing) + skyvar
                spec_box_0 = spectrum.sum(axis = 0)                            # initial box-extracted spectrum
                var_box_0 = var.sum(axis = 0)                                  # initial variance guess
                # Mnew = np.ones_like(M[max(min(peaks) - meta.window, 0):min(max(peaks) + meta.window, rmax),:])
                # Mnew = M[max(peaks_mid - 110, 0):min(peaks_mid + 110, rmax),:]
                Mnew = M[max(min(peaks) - meta.window, 0):min(max(peaks) + meta.window, rmax),:]
                # TODO: Just use meta to reduce the number of parameters which are given to optextr
                if meta.opt_extract==True: [f_opt_0, var_opt_0, numoutliers] = optextr.optextr(spectrum, err, spec_box_0, var_box_0, Mnew, meta.nsmooth, meta.sig_cut, meta.save_optextr_plot, i, ii, meta)
                else: [f_opt, var_opt] = [spec_box_0,var_box_0]

                # NOTE: Sums up spectra and variance for all the differenced images
                spec_opt += f_opt_0
                var_opt += var_opt_0
                spec_box += spec_box_0
                var_box += var_box_0

        ######################################################################################################################################
        # TODO: Q: int(meta.refpix[orbnum, 1]) + meta.LTV1 is kinda sus
        # TODO: Q: in util.get_wave_grid we have:
        # TODO: Q: disp_solution = geo.dispersion(meta.refpix[i,1], -meta.LTV2+j)
        # TODO: Q: delx = 0.5 + np.arange(meta.subarray_size) - (meta.refpix[i,2] + meta.LTV1 + meta.POSTARG1/meta.platescale)

        template_waves = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]
        # print(template_waves)

        # NOTE: Corrects for wavelength drift over time
        if meta.correct_wave_shift == True:
            if i in meta.new_visit_idx_sp:
                # NOTE: Store x & y data if it's the first exposure in the visit
                if meta.correct_wave_shift_refspec == True:
                    x_data_firstexpvisit, y_data_firstexpvisit, leastsq_res = util.correct_wave_shift_fct_0(meta, orbnum, cmin, cmax, spec_opt, i)
                else:
                    x_data_firstexpvisit, y_data_firstexpvisit, leastsq_res = util.correct_wave_shift_fct_00(meta, orbnum, cmin, cmax, spec_opt, i)
                wvls = np.copy(x_data_firstexpvisit)
                if (meta.save_drift_plot or meta.show_drift_plot) and (not meta.correct_wave_shift_refspec):
                    leastsq_res_all.append(leastsq_res)
            else:
                wvls, leastsq_res = util.correct_wave_shift_fct_1(meta, orbnum, cmin, cmax, spec_opt, x_data_firstexpvisit, y_data_firstexpvisit, i)
                if meta.save_drift_plot or meta.show_drift_plot:
                    leastsq_res_all.append(leastsq_res)
        # NOTE: If you dont want to correct it:
        else:
            wvls = template_waves

        # NOTE: Stores 1d spectra into list for plot
        if meta.opt_extract and (meta.save_sp1d_diff_plot or meta.show_sp1d_diff_plot):
            sp1d_all.append(np.interp(wvl_hires, wvls, spec_opt))
        if not meta.opt_extract and (meta.save_sp1d_diff_plot or meta.show_sp1d_diff_plot):
            sp1d_all.append(np.interp(wvl_hires, wvls, spec_box))

        # NOTE: Plot of the 1d spectrum
        if meta.save_sp1d_plot or meta.show_sp1d_plot:
            if meta.opt_extract:
                plots.sp1d(wvls, spec_box, meta, i, spec_opt=spec_opt)
            else:
                plots.sp1d(wvls, spec_box, meta, i)

        # NOTE: Adds rows to the astropy tables
        table_white.add_row([meta.t_mjd_sp[i], meta.t_bjd_sp[i], meta.t_visit_sp[i], meta.t_orbit_sp[i], visnum, orbnum, scan, sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box)])
        n = len(spec_opt)
        for ii in np.arange(n):
            table_spec.add_row([meta.t_mjd_sp[i], meta.t_bjd_sp[i], meta.t_visit_sp[i], meta.t_orbit_sp[i], visnum, orbnum, scan, spec_opt[ii], var_opt[ii], wvls[ii]])
        table_diagnostics.add_row([nspectra, meta.t_mjd_sp[i], numoutliers, skymedian, sum(np.isnan(spec_opt))])
        nspectra += 1
        print('\n')

    # NOTE: Save results in the astropy tables
    if meta.output:
        ascii.write(table_white, dirname / 'lc_white.txt',
                    format='ecsv', overwrite=True)
        ascii.write(table_spec, dirname / 'lc_spec.txt',
                    format='ecsv', overwrite=True)
        ascii.write(table_diagnostics, dirname / 'diagnostics.txt',
                    format='ecsv', overwrite=True)
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir / f'WFC3_{meta.eventlabel}_Meta_Save', save=[])

    # NOTE: Make Plots
    if meta.save_bkg_evo_plot or meta.show_bkg_evo_plot:
        plots.bkg_evo(bkg_evo, meta)

    if meta.save_sp1d_diff_plot or meta.show_sp1d_diff_plot:
        sp1d_all = np.array(sp1d_all)
        sp1d_all_diff = np.diff(sp1d_all, axis=0)
        plots.sp1d_diff(sp1d_all_diff, meta, wvl_hires)

    if meta.save_utr_aper_evo_plot or meta.show_utr_aper_evo_plot:
        plots.utr_aper_evo(peaks_all, meta)

    if meta.save_drift_plot or meta.show_drift_plot:
        plots.drift(leastsq_res_all, meta)

    if meta.save_rowshift_stretch_plot or meta.show_rowshift_stretch_plot:
        plots.rowshift_stretch(meta)

    print('Finished s20 \n')
    return meta
