import os, glob, scipy, numpy

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from numpy import *
import scipy.signal
from pylab import *
from ..lib import optextr
from scipy.optimize import leastsq
from astropy.io import ascii
from datetime import datetime
from astropy.table import QTable
from ..lib import manageevent as me
from ..lib import util
from ..lib import plots
from scipy.signal import find_peaks
from tqdm import tqdm


def run20(eventlabel, workdir, meta=None):

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    # load in more information into meta
    meta = util.ancil(meta, s20=True)

    if meta.grism == 'G102':
        from ..lib import geometry102 as geo
    elif meta.grism == 'G141':
        from ..lib import geometry as geo
    #else:
    #    print('Error: GRISM in obs_par.cf is neither G102 nor G141!')

    # TODO: cf mentioning


    ###############################################################################################################################################################
    ###############################################################################################################################################################

    #STEP 0: User-set parameters and constants

    dirname = meta.workdir + "/extracted_lc/" + datetime.strftime(datetime.now(), '%Y-%m-%d_%H_%M')
    if not os.path.exists(dirname): os.makedirs(dirname)

    table_white = QTable(names=('t_mjd', 't_bjd', 't_visit','t_orbit', 'ivisit', 'iorbit', 'scan', 'spec_opt', 'var_opt', 'spec_box', 'var_box'))
    table_spec = QTable(names=('t_mjd', 't_bjd', 't_visit','t_orbit', 'ivisit', 'iorbit', 'scan', 'spec_opt', 'var_opt', 'template_waves'))
    table_diagnostics = QTable(names=('nspectra', 't_mjd', 'numoutliers', 'skymedian', 'shift', "# nans"))


    # Only use the files which are in the visits one is interested in. These are set by setting which_visits in obs_par
    files_sp = meta.files_sp
    meta.nexp = len(files_sp)
    nspectra = 0                                                #iterator variable to track number of spectra reduced

    peaks_all = []
    bkg_lc = []

    cmin_list= []
    cmax_list= []

    spec_opt_interp_all = []
    wvl_hires = np.linspace(10000, 17800, 1000)

    spec1d_all = []

    print('in total visit, orbit:', (meta.norbit, meta.nvisit), '\n')

    for i in tqdm(np.arange(len(files_sp), dtype=int)):#
        f = files_sp[i]
        #print("\nProgress: {0}/{1}".format(i+1, len(files_sp)))
        print("Filename: {0}".format(f))
        d = fits.open(f)                                        #opens the file
        scan = meta.scans_sp[i]

        # Plot with good visible background
        if meta.save_spectrum2d_plot or meta.show_spectrum2d_plot:
            plots.spectrum2d(d, meta, i)

        visnum, orbnum = meta.ivisit_sp[i], meta.iorbit_sp_com[i]     #current visit and orbit number
        print('current visit, orbit:', (visnum, orbnum))

        # Plot trace
        if meta.save_trace_plot or meta.show_trace_plot:
            plots.plot_trace(d, meta, visnum, orbnum, i)

        offset = meta.offset
        cmin = int(meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale) + meta.BEAMA_i + meta.LTV1 + offset                      #determines left column for extraction (beginning of the trace)
        cmax = min(int(meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale) + meta.BEAMA_f + meta.LTV1 - offset, meta.subarray_size)     #right column (end of trace, or edge of detector)
        rmin, rmax = int(meta.rmin), int(meta.rmax)                     #top and bottom row for extraction (specified in obs_par.txt)
        #print(cmin, cmax, rmin, rmax)
        #skycmin, skycmax, skyrmin, skyrmax = meta.skycmin, meta.skycmax, meta.skyrmin, meta.skyrmax

        #D = np.zeros_like(d[1].data[rmin:rmax,cmin:cmax])                        #array to store the background-subtracted data
        #outlier_array = np.zeros_like(d[1].data[rmin:rmax,cmin:cmax])                    #array used to determine which pixel is the biggest outlier
        M = np.ones_like(d[1].data[rmin:rmax, cmin:cmax])                        #mask for bad pixels

        meta.wave_grid = util.get_wave_grid(meta)  # gets grid of wavelength solutions for each orbit and row
        flatfield = util.get_flatfield(meta)

        bpix = d[3].data[rmin:rmax,cmin:cmax]
        badpixind =  (bpix==4)|(bpix==512)|(flatfield[orbnum][rmin:rmax, cmin:cmax] == -1.)    #selects bad pixels
        M[badpixind] = 0.0                                        #initializes bad pixel mask

        spec_box = np.zeros(cmax - cmin)                                #box extracted standard spectrum
        spec_opt = np.zeros(cmax - cmin)                                #optimally extracted spectrum
        var_box = np.zeros(cmax - cmin)                                #box spectrum variance
        var_opt = np.zeros(cmax - cmin)                                #optimal spectrum variance


        #########################################################################################################################################################
        # loops over up-the-ramp-samples (skipping first two very short exposures); gets all needed input for optextr routine                    #
        #########################################################################################################################################################

        for ii in range(d[0].header['nsamp']-2):
            print("Progress (up-the-ramp-samples): {0}/{1}".format(ii+1, d[0].header['nsamp']-2))
            diff = d[ii*5 + 1].data[rmin:rmax,cmin:cmax] - d[ii*5 + 6].data[rmin:rmax,cmin:cmax]    #creates image that is the difference between successive scans

            #diff = diff/flatfield[orbnum][rmin:rmax, cmin:cmax]                               #flatfields the differenced image

            #idx = np.argmax(scipy.signal.medfilt(np.sum(diff, axis = 1),3))                         #computes spatial index of peak counts in the difference image

            rowsum = np.sum(diff, axis=1)                   # sum of every row
            rowsum_absder = abs(rowsum[1:] - rowsum[:-1])   # absolute derivative

            peaks, _ = find_peaks(rowsum_absder, height=max(rowsum_absder * 0.2), distance=meta.window)
            print('peak positions: ', peaks)
            peaks = peaks[:2]
            #idx = int(np.mean(peaks))

            peaks_all.append(peaks)

            #estimates sky background and variance
            fullframe_diff = d[ii*5 + 1].data - d[ii*5 + 6].data                                       #fullframe difference between successive scans
            #fullframe_diff = d[1].data

            #if meta.background_box:
            #    skymedian = np.median(fullframe_diff[skyrmin:skyrmax,skycmin:skycmax])                    #estimates the background counts
            #    skyvar = util.median_abs_dev(fullframe_diff[skyrmin:skyrmax,skycmin:skycmax].flatten())            #variance for the background count estimate

            below_threshold = fullframe_diff<meta.background_thld
            skymedian = np.median(fullframe_diff[below_threshold].flatten())  # estimates the background counts
            bkg_lc.append(skymedian)
            skyvar = util.median_abs_dev(fullframe_diff[below_threshold].flatten())  # variance for the background count estimate
            if meta.save_bkg_hist_plot or meta.show_bkg_hist_plot:
                plots.bkg_hist(fullframe_diff, skymedian, meta, i, ii)
            print('bkg: ', skymedian, skyvar)
                #skymedian, skyvar = bkg_histogram_median

            diff = diff - skymedian                                    #subtracts the background

            #interpolation to correct for bad pixels and the fact that the wavelength solution changes row by row
            """for jj in range(rmax):
                goodidx = M[jj,:] == 1.0                            #selects good pixels
                diff[jj,:] = np.interp(ancil.wave_grid[0, int(ancil.refpix[0,1]) + ancil.LTV1, \
                    cmin:cmax], ancil.wave_grid[0, jj, cmin:cmax][goodidx], diff[jj,goodidx])    #LK interpolation 8/18"""


            spectrum = diff[max(min(peaks) - meta.window, 0):min(max(peaks) + meta.window, rmax),:]        #selects postage stamp centered around spectrum

            if meta.save_uptheramp_plot or meta.show_uptheramp_plot:
                plots.uptheramp(diff, meta, i, ii, orbnum, rowsum, rowsum_absder, peaks)


            #stores median of column that has spectrum on it (+/- 5 pix from center) for ii = 0

            err = np.zeros_like(spectrum) + float(meta.rdnoise)**2 + skyvar
            var = abs(spectrum) + float(meta.rdnoise)**2 +skyvar                #variance estimate: Poisson noise from photon counts (first term)  + readnoise (factor of 2 for differencing) + skyvar
            spec_box_0 = spectrum.sum(axis = 0)                            #initial box-extracted spectrum
            var_box_0 = var.sum(axis = 0)                                #initial variance guess
            newM = np.ones_like(spectrum)  # not masking any pixels because we interpolated over them
            #print(M.shape)
            #print(newM.shape)
            if meta.opt_extract==True: [f_opt_0, var_opt_0, numoutliers] = optextr.optextr(spectrum, err, spec_box_0, var_box_0, newM, meta.nsmooth, meta.sig_cut, meta.diagnostics)
            else: [f_opt, var_opt] = [spec_box_0,var_box_0]

            #sums up spectra and variance for all the differenced images
            spec_opt += f_opt_0
            var_opt += var_opt_0
            spec_box += spec_box_0
            var_box += var_box_0

        ######################################################################################################################################
        #print(d[0].header['expstart'])
        #time = (d[0].header['expstart'] + d[0].header['expend'])/(2.0) + 2400000.5                    #converts time to BJD_TDB; see Eastman et al. 2010 equation 4
        #time = time + (32.184)/(24.0*60.0*60.0)
        #print("visnum, suntimecorr (minutes)", visnum, (suntimecorr.suntimecorr(ancil.ra, ancil.dec, array([time]), ancil.coordtable[visnum], verbose=False))/(60.0))
        #time = [time] + (suntimecorr.suntimecorr(ancil.ra, ancil.dec, array([time]), ancil.coordtable[visnum], verbose=False))/(60.0*60.0*24.0)

        #phase = (time-ancil.t0)/ancil.period - math.floor((time-ancil.t0)/ancil.period)
        #if phase > 0.5: phase = phase - 1.0                                    #calculates orbital phase

        print(len(np.arange(cmin, cmax)))
        template_waves = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]
        print(len(template_waves))

        np.savetxt(meta.workdir + 'template_waves.txt', list(zip(np.arange(cmin, cmax), template_waves)))

        phase = 0
        shift = 0.
        #corrects for wavelength drift over time
        if meta.correct_wave_shift == True:
            if i in meta.new_visit_idx_sp:
                #if nspectra == 0:
                template_waves = meta.wave_grid[0, int(meta.refpix[orbnum,1]) + meta.LTV1, cmin:cmax]             #LK interpolation 8/18 #use stellar model instead
                g102mask = template_waves > 8200
                #print(template_waves[-1]-template_waves[0])
                #print(cmax-cmin)
                #print((template_waves[-1] - template_waves[0])/(cmax-cmin))
                #https: // matthew - brett.github.io / teaching / smoothing_intro.html
                refspec = np.loadtxt(meta.workdir + '/ancil/bandpass/refspec_v{0}.txt'.format(visnum)).T
                x_refspec, y_refspec_raw = refspec[0]*10000, refspec[1]/max(refspec[1])

                print(x_refspec[-1] - x_refspec[0])

                sigma = 20#1*46.17#0.004*10000
                y_refspec_kernel = np.zeros(y_refspec_raw.shape)

                for x_position, x_position_val in enumerate(x_refspec):
                    kernel = np.exp(-(x_refspec - x_position_val) ** 2 / (2 * sigma ** 2))
                    kernel = kernel / sum(kernel)
                    y_refspec_kernel[x_position] = sum(y_refspec_raw * kernel)

                y_refspec_kernel = y_refspec_kernel/max(y_refspec_kernel)

                #FIXME SZ make this nicer
                x_refspec_new = np.concatenate((np.linspace(-10000, min(x_refspec), 100, endpoint=False),
                                         x_refspec,
                                         np.linspace(max(x_refspec) + 10, 100000, 100, endpoint=False)))
                y_refspec_kernel_new = np.concatenate((np.zeros(100),
                                         y_refspec_kernel,
                                         np.zeros(100)))

                x_data = template_waves[g102mask]
                y_data = (spec_opt/max(spec_opt))[g102mask]

                p0 = [0, 1, 1]
                leastsq_res = leastsq(util.residuals2, p0, args=(x_refspec_new, y_refspec_kernel_new, x_data, y_data))[0]
                print('leastsq_res', leastsq_res)

                if meta.save_refspec_comp_plot or meta.show_refspec_comp_plot:
                    plots.refspec_comp(x_refspec, y_refspec_raw, x_refspec_new, y_refspec_kernel_new, p0, x_data, y_data, leastsq_res, meta, i)

                template_waves = leastsq_res[0] + template_waves * leastsq_res[1]

                #for all other but first exposure in visit exposures
                template_waves_ref = np.copy(template_waves)
                y_data_firstexpvisit = np.copy(spec_opt/max(spec_opt))

                #template_waves -= 70.               #LK hack to get the wavelength solution right
                #print("LK adjusting wavelength solution BY HAND!!!!")
                #shift = 0.
                #template_spectrum = spec_opt                                #makes the first exposure the template spectrum #stellar model * grism throughput
                #best_spec = spec_opt
                #best_var = var_opt

                # else:
                #     #shifts spectrum so it matches the template
                #     [best_spec, best_var, shift] = interpolate_spectrum(spec_opt, np.sqrt(var_opt), template_spectrum, template_waves)
                #     best_var = best_var**2
                #
                #     spec_opt = best_spec                                    #saves the interpolated spectrum
                #     var_opt = best_var
            else:

                #FIXME SZ make this nicer
                x_model = np.concatenate((np.linspace(-10000, min(template_waves_ref), 100, endpoint=False),
                                         template_waves_ref,
                                         np.linspace(max(template_waves_ref) + 10, 100000, 100, endpoint=False)))
                y_model = np.concatenate((np.zeros(100),
                                         y_data_firstexpvisit,
                                         np.zeros(100)))

                x_data = meta.wave_grid[0, int(meta.refpix[orbnum,1]) + meta.LTV1, cmin:cmax]
                y_data = spec_opt/max(spec_opt)

                p0 = [0, 1, 1]
                leastsq_res = leastsq(util.residuals2, p0, args=(x_model, y_model, x_data, y_data))[0]
                print('leastsq_res', leastsq_res)

               # if meta.save_refspec_comp_plot or meta.show_refspec_comp_plot:
               #     plots.refspec_comp(x_vals, y_vals, modelx, modely, p0, datax, datay, leastsq_res, meta, i)
                plots.refspec_comp2(x_model, y_model, p0, x_data, y_data, leastsq_res, meta, i)

                template_waves = leastsq_res[0] + x_data * leastsq_res[1]

        else:
            template_waves = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]


        cmin_list.append(cmin)
        cmax_list.append(cmax)

        spec_opt_interp_all.append(np.interp(wvl_hires, template_waves, spec_opt))

        #wavelengthsolutionoffset = 105.

        #plots wavelength solution compared to a stellar template
        """plt.plot(template_waves - wavelengthsolutionoffset, best_spec/np.max(best_spec), label = 'LK spectrum')
        np.save("temp_spectrum", [template_waves, best_spec])
    
        synth = np.load("syntehtic_spec.npz.npy")
        plt.plot(synth[0], 1.05*synth[1]/np.max(synth[1]), label = 'template')
        plt.legend()
        plt.show()"""

        if meta.save_spectrum1d_spec_opt_plot or meta.show_spectrum1d_spec_opt_plot:
            plots.spectrum1d_spec_opt(cmin,cmax,template_waves, spec_opt, meta, i)

        if meta.save_spectrum1d_spec_box_plot or meta.show_spectrum1d_spec_box_plot:
            plots.spectrum1d_spec_box(cmin,cmax,template_waves, spec_box, meta, i)

        #print(phase, sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time, visnum, orbnum, scan)

       # print(phase[0], sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time[0], visnum, orbnum, scan)
        if meta.output == True:
            table_white.add_row([meta.t_mjd_sp[i], meta.t_bjd_sp[i], meta.t_visit_sp[i], meta.t_orbit_sp[i], visnum, orbnum, scan, sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box)])
            n = len(spec_opt)
            for ii in arange(n):
                table_spec.add_row([meta.t_mjd_sp[i], meta.t_bjd_sp[i], meta.t_visit_sp[i], meta.t_orbit_sp[i], visnum, orbnum, scan, spec_opt[ii], var_opt[ii], template_waves[ii]])
            #print(nspectra, time[0], numoutliers, skymedian, shift, file=diagnosticsfile)
            table_diagnostics.add_row([nspectra, meta.t_mjd_sp[i], numoutliers, skymedian, shift, sum(np.isnan(spec_opt))])

        nspectra += 1
        #if nspectra%10 == 0: print("Extraction", '{0:1f}'.format(float(nspectra)/float(len(ancil.files))*100.), "% complete, time elapsed (min) =", '{0:0.1f}'.format(clock()/60.))




        spec1d_all.append(spec_opt)


        #print("sky background!!", skymedian, visnum)
        #print("length of spectrum", len(spec_opt))
        print("# nans", sum(np.isnan(spec_opt)))
        print(nspectra, meta.t_bjd_sp[i], sum(spec_opt), np.sum(spec_box), visnum, f, shift)
        print('\n')

    spec1d_all_diff = np.diff(spec1d_all, axis=0)
    for iiii in range(len(spec1d_all_diff)):
        plt.plot(range(len(spec1d_all_diff[iiii])), spec1d_all_diff[iiii])
        if not os.path.isdir(meta.workdir + '/figs/scale/'):
            os.makedirs(meta.workdir + '/figs/scale/')
        plt.savefig(meta.workdir + '/figs/scale/scale{0}.png'.format(iiii))
        plt.close()

    if meta.output == True:
        ascii.write(table_white, dirname+'/lc_white.txt', format='ecsv', overwrite=True)
        ascii.write(table_spec, dirname+'/lc_spec.txt', format='ecsv', overwrite=True)
        ascii.write(table_diagnostics, dirname+'/diagnostics.txt', format='ecsv', overwrite=True)
    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    if meta.save_spectrum1d_spec_opt_diff_plot or meta.show_spectrum1d_spec_opt_diff_plot:
        spec_opt_interp_all = np.array(spec_opt_interp_all)
        spec_opt_interp_all_diff = np.diff(spec_opt_interp_all, axis=0)
        plots.spectrum1d_spec_opt_diff_plot(spec_opt_interp_all_diff, meta, wvl_hires)

    if meta.save_uptheramp_plot_total or meta.show_uptheramp_plot_total:
        plots.bkg_lc(bkg_lc, meta)
        plots.uptheramp_evolution(peaks_all, meta)

    if meta.save_c_diag_plot or meta.show_c_diag_plot:
        plots.c_diag(cmin_list, cmax_list, meta)

    print("I made a change to how the wavelength interpolation is being done")
    print("it's commented as LK 8/18")
    print("Should check at some point to see if it's right")

    print('Finished s20 \n')

    return meta

    # - create new directory
    # - create txt files in there to save data

    # - start loop over spectra
    # - save a picture of spectrum
    # - save a plot of the trace
    # - rmin:rmax,cmin:cmax:
    # -- determines left column for extraction (beginning of the trace)
    # -- right column (end of trace, or edge of detector)
    # -- top and bottom row for extraction (specified in obs_par.txt)
    # - mask for bad pixels
    # - save picture of spectrum with extraction box and background box
    # - new loop over scans and take differences
    # - diff/flatfield
    # - computes spatial index of peak counts in the difference image
    # - calculate median background in fullframe difference image
    # - diff = diff - skymedian
    # - optimal extraction
    # - sums up spectra and variance for all the differenced images
    # - end of loop over scans
    # - #corrects for wavelength drift over time using the template
    # - save data
