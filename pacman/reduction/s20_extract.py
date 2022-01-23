import os, glob, scipy, numpy, sys
import numpy as np
from astropy.io import ascii, fits
#from numpy import *
#from pylab import *
from ..lib import optextr
from scipy.optimize import leastsq
from datetime import datetime
from astropy.table import QTable
from scipy.signal import find_peaks
from tqdm import tqdm
from ..lib import manageevent as me
from ..lib import util
from ..lib import plots

def run20(eventlabel, workdir, meta=None):
    """
    This function extracts the spectrum and saves the total flux and the flux as a function of wavelength into files.
    """

    print('Starting s20')

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    # load in more information into meta
    meta = util.ancil(meta, s20=True)

    ###############################################################################################################################################################
    #STEP 0: Set up files and directories
    ###############################################################################################################################################################

    dirname = meta.workdir + "/extracted_lc/" + datetime.strftime(datetime.now(), '%Y-%m-%d_%H-%M-%S')
    if not os.path.exists(dirname): os.makedirs(dirname)

    # initialize the astropy tables where we will save the extracted spectra
    if meta.output == True:
        table_white = QTable(names=('t_mjd', 't_bjd', 't_visit','t_orbit', 'ivisit', 'iorbit', 'scan', 'spec_opt', 'var_opt', 'spec_box', 'var_box'))
        table_spec = QTable(names=('t_mjd', 't_bjd', 't_visit','t_orbit', 'ivisit', 'iorbit', 'scan', 'spec_opt', 'var_opt', 'template_waves'))
        table_diagnostics = QTable(names=('nspectra', 't_mjd', 'numoutliers', 'skymedian', 'shift', "# nans"))

    files_sp = meta.files_sp     # spectra files
    meta.nexp = len(files_sp)    # number of exposures
    nspectra = 0                 # iterator variable to track number of spectra reduced

    # the following lists are used for diagnostic plots
    if meta.save_utr_aper_evo_plot or meta.show_utr_aper_evo_plot:
        peaks_all = []
    if meta.save_bkg_evo_plot or meta.show_bkg_evo_plot:
        bkg_evo = []
    if meta.save_sp1d_diff_plot or meta.show_sp1d_diff_plot:
        sp1d_all = []
        wvl_hires = np.linspace(7000, 17800, 1000)




    print('in total #visits, #orbits:', (meta.nvisit, meta.norbit), '\n')
    # in order to have the correct order of print() with tqdm, i added file=sys.stdout
    # source: https://stackoverflow.com/questions/36986929/redirect-print-command-in-python-script-through-tqdm-write
    for i in tqdm(np.arange(2, dtype=int), desc='***************** Looping over files', file=sys.stdout):#tqdm(np.arange(len(files_sp), dtype=int)):
        f = files_sp[i]                     # current file
        print("\nFilename: {0}".format(f))
        d = fits.open(f)                    # opens the file
        scan = meta.scans_sp[i]             # scan direction of the spectrum.

        # Plot with good visible background
        if meta.save_sp2d_plot or meta.show_sp2d_plot:
            plots.sp2d(d, meta, i)

        visnum, orbnum = meta.ivisit_sp[i], meta.iorbit_sp_cumulative[i]     #current visit and cumulative orbit number
        print('current visit, orbit: ', (visnum, orbnum))

        # Plot trace
        #TODO: Q: What is the meaning of the y location of the trace?
        if meta.save_trace_plot or meta.show_trace_plot:
            plots.trace(d, meta, visnum, orbnum, i)

        #TODO: Q: What is the offset used for?
        offset = meta.offset
        #TODO: SPEED UP: calculation of the start and end of the trace could be moved to util.py. It's also used in plots.plot_trace. Also in plots.utr
        #TODO: Q: There's also cmin and cmax in the pcf file. Keep?
        cmin = int(meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale) + meta.BEAMA_i + meta.LTV1 + offset                      #determines left column for extraction (beginning of the trace)
        cmax = min(int(meta.refpix[orbnum,2] + meta.POSTARG1/meta.platescale) + meta.BEAMA_f + meta.LTV1 - offset, meta.subarray_size)     #right column (end of trace, or edge of detector)
        #TODO: Q: Should we leave the decision on the top and bottom row to the user or just use minimum and max possible rows
        #TODO: Q: Handbook: First-order spectra for both the G102 and G141 grisms comfortably fit within 512 × 512 and 256 × 256 pixel subarrays.
        #TODO: Q: Why are the dimensions 266x266?
        rmin, rmax = int(meta.rmin), int(meta.rmax)                     #top and bottom row for extraction (specified in obs_par.txt)

        #TODO: Q: Should we keep a way for the user to decide on a background box?
        #skycmin, skycmax, skyrmin, skyrmax = meta.skycmin, meta.skycmax, meta.skyrmin, meta.skyrmax

        #D = np.zeros_like(d[1].data[rmin:rmax,cmin:cmax])                        #array to store the background-subtracted data
        #outlier_array = np.zeros_like(d[1].data[rmin:rmax,cmin:cmax])                    #array used to determine which pixel is the biggest outlier


        #TODO: SPEED UP: This calculation is done for every file again from new!
        #TODO: Move it to util.py so its done just once
        meta.wave_grid = util.get_wave_grid(meta)  # gets grid of wavelength solutions for each orbit and row

        #TODO: Q: The flat field file is used to mask bad pixels.
        #TODO: Q: "sets flatfield pixels below 0.5 to -1 so they can be masked"
        #TODO: Q: Keep that?

        #TODO: Q: M is not used
        #TODO: Q: Earlier M was used to do the interpolation to correct for bad pixels (line 154)
        M = np.ones_like(d[1].data[rmin:rmax, cmin:cmax])                        #mask for bad pixels
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
        # in order to not print a new line with tqdm every time, I added leave=True, position=0
        # source: https://stackoverflow.com/questions/41707229/tqdm-printing-to-newline
        for ii in tqdm(np.arange(d[0].header['nsamp']-2, dtype=int), desc='--- Looping over up-the-ramp-samples', leave=True, position=0):
            #print("Progress (): {0}/{1}".format(ii+1, d[0].header['nsamp']-2))
            diff = d[ii*5 + 1].data[rmin:rmax,cmin:cmax] - d[ii*5 + 6].data[rmin:rmax,cmin:cmax]    #creates image that is the difference between successive scans

            #diff = diff/flatfield[orbnum][rmin:rmax, cmin:cmax]                               #flatfields the differenced image

            # determine the aperture cutout
            rowsum = np.sum(diff, axis=1)                   # sum of every row
            rowsum_absder = abs(rowsum[1:] - rowsum[:-1])   # absolute derivative in order to determine where the spectrum is
            #we use scipy.signal.find_peaks to determine in which rows the spectrum starts and ends
            #TODO: Think about better values for height and distance
            #TODO: Finding the row with the highest change in flux (compared to row above and below) isnt robust against outliers!
            #TODO: Maybe instead use the median flux in a row??
            peaks, _ = find_peaks(rowsum_absder, height=max(rowsum_absder * 0.2), distance=5)
            #print('peak positions: ', peaks)
            peaks = peaks[:2] # only take the two biggest peaks (there should be only 2)

            #stores the locations of the peaks for every file and up-the-ramp-samples
            if meta.save_utr_aper_evo_plot or meta.show_utr_aper_evo_plot:
                peaks_all.append(peaks)

            #estimates sky background and variance
            fullframe_diff = d[ii*5 + 1].data - d[ii*5 + 6].data                                       #fullframe difference between successive scans

            #if meta.background_box:
            #    skymedian = np.median(fullframe_diff[skyrmin:skyrmax,skycmin:skycmax])                    #estimates the background counts
            #    skyvar = util.median_abs_dev(fullframe_diff[skyrmin:skyrmax,skycmin:skycmax].flatten())            #variance for the background count estimate

            ### BACKGROUND SUBTRACTION
            below_threshold = fullframe_diff < meta.background_thld # mask with all pixels below the user defined threshold
            skymedian = np.median(fullframe_diff[below_threshold].flatten())  # estimates the background counts by taking the flux median of the pixels below the flux threshold
            if meta.save_bkg_evo_plot or meta.show_bkg_evo_plot:
                bkg_evo.append(skymedian)
            skyvar = util.median_abs_dev(fullframe_diff[below_threshold].flatten())  # variance for the background count estimate
            if meta.save_bkg_hist_plot or meta.show_bkg_hist_plot:
                plots.bkg_hist(fullframe_diff, skymedian, meta, i, ii)
            diff = diff - skymedian                                    #subtracts the background

            #TODO: Q: ?
            #interpolation to correct for bad pixels and the fact that the wavelength solution changes row by row
            """for jj in range(rmax):
                goodidx = M[jj,:] == 1.0                            #selects good pixels
                diff[jj,:] = np.interp(ancil.wave_grid[0, int(ancil.refpix[0,1]) + ancil.LTV1, \
                    cmin:cmax], ancil.wave_grid[0, jj, cmin:cmax][goodidx], diff[jj,goodidx])    #LK interpolation 8/18"""

            # selects postage stamp centered around spectrum
            # we use a bit more data by using the user defined window
            spectrum = diff[max(min(peaks) - meta.window, 0):min(max(peaks) + meta.window, rmax),:]

            if meta.save_utr_plot or meta.show_utr_plot:
                plots.utr(diff, meta, i, ii, orbnum, rowsum, rowsum_absder, peaks)

            err = np.zeros_like(spectrum) + float(meta.rdnoise)**2 + skyvar
            var = abs(spectrum) + float(meta.rdnoise)**2 + skyvar                #variance estimate: Poisson noise from photon counts (first term)  + readnoise (factor of 2 for differencing) + skyvar
            spec_box_0 = spectrum.sum(axis = 0)                            #initial box-extracted spectrum
            var_box_0 = var.sum(axis = 0)                                #initial variance guess
            newM = np.ones_like(spectrum)  # not masking any pixels because we interpolated over them
            #print(M.shape)
            #print(newM.shape)
            #TODO: Just use meta to reduce the number of parameters which are given to optextr
            if meta.opt_extract==True: [f_opt_0, var_opt_0, numoutliers] = optextr.optextr(spectrum, err, spec_box_0, var_box_0, newM, meta.nsmooth, meta.sig_cut, meta.save_optextr_plot, i, ii, meta)
            else: [f_opt, var_opt] = [spec_box_0,var_box_0]

            #sums up spectra and variance for all the differenced images
            spec_opt += f_opt_0
            var_opt += var_opt_0
            spec_box += spec_box_0
            var_box += var_box_0

        ######################################################################################################################################
        #TODO: Q: int(meta.refpix[orbnum, 1]) + meta.LTV1 is kinda sus
        #TODO: Q: in util.get_wave_grid we have:
        #TODO: Q: disp_solution = geo.dispersion(meta.refpix[i,1], -meta.LTV2+j)
        #TODO: Q: delx = 0.5 + np.arange(meta.subarray_size) - (meta.refpix[i,2] + meta.LTV1 + meta.POSTARG1/meta.platescale)

        template_waves = meta.wave_grid[0, int(meta.refpix[orbnum, 1]) + meta.LTV1, cmin:cmax]

        #TODO: Q: shift removable right?
        shift = 0.

        #corrects for wavelength drift over time
        if meta.correct_wave_shift == True:
            if i in meta.new_visit_idx_sp:
                #if nspectra == 0:
                template_waves = meta.wave_grid[0, int(meta.refpix[orbnum,1]) + meta.LTV1, cmin:cmax]             #LK interpolation 8/18 #use stellar model instead
                g102mask = template_waves > 8200 # we dont use the spectrum below 8200 angstrom for the interpolation as the reference bandpass cuts out below this wavelength

                #https: // matthew - brett.github.io / teaching / smoothing_intro.html
                refspec = np.loadtxt(meta.workdir + '/ancil/refspec/refspec.txt').T
                x_refspec, y_refspec_raw = refspec[0]*1e10, refspec[1]/max(refspec[1])

                kernel_bool = False
                if kernel_bool:
                    print(x_refspec[-1] - x_refspec[0])

                    sigma = 20#1*46.17#0.004*10000
                    y_refspec_kernel = np.zeros(y_refspec_raw.shape)

                    for x_position, x_position_val in enumerate(x_refspec):
                        kernel = np.exp(-(x_refspec - x_position_val) ** 2 / (2 * sigma ** 2))
                        kernel = kernel / sum(kernel)
                        y_refspec_kernel[x_position] = sum(y_refspec_raw * kernel)
                    y_refspec_kernel = y_refspec_kernel/max(y_refspec_kernel)

                y_refspec_kernel = y_refspec_raw
                #FIXME SZ make this nicer
                x_refspec_new = np.concatenate((np.linspace(-10000, min(x_refspec), 100, endpoint=False),
                                         x_refspec,
                                         np.linspace(max(x_refspec) + 10, 100000, 100, endpoint=False)))
                y_refspec_kernel_new = np.concatenate((np.zeros(100),
                                         y_refspec_kernel,
                                         np.zeros(100)))

                x_data = template_waves[g102mask]
                y_data = (spec_opt/max(spec_opt))[g102mask]

                #np.savetxt(meta.workdir + '/first_exp.txt', list(zip(x_data, y_data)))

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

        # stores 1d spectra into list for plot
        if meta.opt_extract and (meta.save_sp1d_diff_plot or meta.show_sp1d_diff_plot):
            sp1d_all.append(np.interp(wvl_hires, template_waves, spec_opt))
        if not meta.opt_extract and (meta.save_sp1d_diff_plot or meta.show_sp1d_diff_plot):
            sp1d_all.append(np.interp(wvl_hires, template_waves, spec_box))

        #wavelengthsolutionoffset = 105.

        #plots wavelength solution compared to a stellar template
        """plt.plot(template_waves - wavelengthsolutionoffset, best_spec/np.max(best_spec), label = 'LK spectrum')
        np.save("temp_spectrum", [template_waves, best_spec])
    
        synth = np.load("syntehtic_spec.npz.npy")
        plt.plot(synth[0], 1.05*synth[1]/np.max(synth[1]), label = 'template')
        plt.legend()
        plt.show()"""

        # plot of the 1d spectrum
        if meta.save_sp1d_plot or meta.show_sp1d_plot:
            if meta.opt_extract:
                plots.sp1d(template_waves, spec_box, meta, i, spec_opt = spec_opt)
            else:
                plots.sp1d(template_waves, spec_box, meta, i)


        #print(phase, sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time, visnum, orbnum, scan)

       # print(phase[0], sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time[0], visnum, orbnum, scan)
        # TODO: Q: Keep asking if user wants output?
        if meta.output == True:
            table_white.add_row([meta.t_mjd_sp[i], meta.t_bjd_sp[i], meta.t_visit_sp[i], meta.t_orbit_sp[i], visnum, orbnum, scan, sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box)])
            n = len(spec_opt)
            for ii in np.arange(n):
                table_spec.add_row([meta.t_mjd_sp[i], meta.t_bjd_sp[i], meta.t_visit_sp[i], meta.t_orbit_sp[i], visnum, orbnum, scan, spec_opt[ii], var_opt[ii], template_waves[ii]])
            #print(nspectra, time[0], numoutliers, skymedian, shift, file=diagnosticsfile)
            table_diagnostics.add_row([nspectra, meta.t_mjd_sp[i], numoutliers, skymedian, shift, sum(np.isnan(spec_opt))])

        nspectra += 1
        #if nspectra%10 == 0: print("Extraction", '{0:1f}'.format(float(nspectra)/float(len(ancil.files))*100.), "% complete, time elapsed (min) =", '{0:0.1f}'.format(clock()/60.))



        #print("sky background!!", skymedian, visnum)
        #print("length of spectrum", len(spec_opt))
        #TODO: Will break if user doesn't use optimal extraction!!
        print("# nans", sum(np.isnan(spec_opt)))
        #print(nspectra, meta.t_bjd_sp[i], sum(spec_opt), np.sum(spec_box), visnum, f, shift)
        print('\n')


    # Save results
    if meta.output == True:
        ascii.write(table_white, dirname+'/lc_white.txt', format='ecsv', overwrite=True)
        ascii.write(table_spec, dirname+'/lc_spec.txt', format='ecsv', overwrite=True)
        ascii.write(table_diagnostics, dirname+'/diagnostics.txt', format='ecsv', overwrite=True)
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])


    # Make Plots
    if meta.save_bkg_evo_plot or meta.show_bkg_evo_plot:
        plots.bkg_evo(bkg_evo, meta)

    if meta.save_sp1d_diff_plot or meta.show_sp1d_diff_plot:
        sp1d_all = np.array(sp1d_all)
        sp1d_all_diff = np.diff(sp1d_all, axis=0)
        plots.sp1d_diff(sp1d_all_diff, meta, wvl_hires)

    if meta.save_utr_aper_evo_plot or meta.show_utr_aper_evo_plot:
        plots.utr_aper_evo(peaks_all, meta)


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
