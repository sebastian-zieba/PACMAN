import sys
sys.path.insert(0, './util')
import os, glob, scipy, numpy
from astropy.io import fits
from numpy import *
import scipy.signal
import matplotlib.pyplot as plt
from pylab import *
import numpy.ma as ma
import suntimecorr
import optextr 
import copy
from scipy.optimize import leastsq
from astropy.io import ascii
from time import clock
from datetime import datetime
from shutil import copyfile
import pickle
import math
import yaml
import ancil
from astropy.table import QTable

from astropy.io import fits

#import obs_par and fit_par
obs_par_path = "config/obs_par.yaml"
with open(obs_par_path, 'r') as file:
    obs_par = yaml.safe_load(file)
fit_par = "config/fit_par.txt"
ancil = ancil.AncillaryData(obs_par, fit_par=fit_par)

if ancil.grism == 'G102':
    import geometry102 as geo
elif ancil.grism == 'G141':
    import geometry as geo
else:
    print('Error: GRISM in obs_par.yaml is neither G102 nor G141!')


### FUNCTIONS

def get_wave_grid(ancil):
    wave_grid = np.empty((ancil.norb*ancil.nvisit, ancil.subarray_size, ancil.subarray_size))

    #calculates wavelength solution row by row for each orbit
    for i in range(ancil.norb*ancil.nvisit):
        for j in range(ancil.subarray_size):
            disp_solution = geo.dispersion(ancil.refpix[i,1], -ancil.LTV2+j)
            delx = 0.5 + np.arange(ancil.subarray_size) - (ancil.refpix[i,2] + ancil.LTV1 + ancil.POSTARG1/ancil.platescale)
            wave_grid[i, j, :] = disp_solution[0] + delx*disp_solution[1]

    return wave_grid

def residuals(params, template_waves, template, spectrum, error):
    shift, scale = params
    fit = scale*np.interp(template_waves, template_waves-shift, spectrum)
    x = (template-fit)/error
    return (template-fit)/error

def interpolate_spectrum(spectrum, error, template, template_waves):
    p0 = [1., 1.0]                                        #initial guess for parameters shift and scale
    plsq, success  = leastsq(residuals, p0, args=(template_waves, template, spectrum, error))
    shift, scale = plsq
    interp_spectrum = np.interp(template_waves, template_waves-shift, spectrum)
    interp_error = np.interp(template_waves, template_waves-shift, error)
    return [interp_spectrum, interp_error, shift]

# def read_dict(filename):
#     #with open(filename,"r") as text:
#     #    return dict(line.strip().split() for line in text)
#     dict = {}
#     with open(filename,"r") as text:
#         for line in text:
#             key, value = line.split()
#             value = str_to_num(value)
#             if value == 'True': value = True
#             if value == 'False': value = False
#             dict[key] = value
#     return dict

def median_abs_dev(vec):
    med = ma.median(vec)
    return ma.median(abs(vec - med))

#def getphase(t):
#    phase = (t - t0)/period
#    return phase - int(phase)

def get_flatfield(ancil):                    #function that flatfields a data array D, which starts at [minr, cmin] of hdu[1].data
    flat = fits.open(ancil.flat)                #reads in flatfield cube
    WMIN = flat[0].header['WMIN']                #constants for computing flatfield coefficients
    WMAX = flat[0].header['WMAX']

    a0 = flat[0].data[-ancil.LTV1:-ancil.LTV1+ancil.subarray_size, -ancil.LTV2:-ancil.LTV2+ancil.subarray_size]
    a1 = flat[1].data[-ancil.LTV1:-ancil.LTV1+ancil.subarray_size, -ancil.LTV2:-ancil.LTV2+ancil.subarray_size]
    a2 = flat[2].data[-ancil.LTV1:-ancil.LTV1+ancil.subarray_size, -ancil.LTV2:-ancil.LTV2+ancil.subarray_size]

    flatfield = []
    for i in range(ancil.norb*ancil.nvisit):
        wave = ancil.wave_grid[i, :]
        x = (wave - WMIN)/(WMAX-WMIN)
        flatfield.append(a0+a1*x+a2*x**2)
        flatfield[i][flatfield[i] < 0.5] = -1.        #sets flatfield pixels below 0.5 to -1 so they can be masked
    return flatfield

def plot_trace(d, ancil):
    trace = geo.trace(ancil.refpix[:,1], ancil.refpix[:,2])                #determines trace coefficients
    #orbnum=ancil.orbnum
    trace_i = ancil.refpix[orbnum,2] + ancil.POSTARG1/ancil.platescale + float(ancil.BEAMA_i + ancil.LTV1)          #start of trace
    trace_f = ancil.refpix[orbnum,2] + ancil.POSTARG1/ancil.platescale + float(ancil.BEAMA_f + ancil.LTV1)          #end of trace
    tracex = np.linspace(trace_i, trace_f,100)                    #x values over which to compute trace
    tracey = ancil.refpix[orbnum,1] + ancil.LTV2 + trace[0][orbnum] + \
    trace[1][orbnum]*(tracex - tracex[0])                    #y values of trace
    plt.imshow(d[1].data, origin = 'lower', vmin=0, vmax=40000)            #plots raw image
    plt.colorbar()
    plt.plot(tracex, tracey, color='yellow', linewidth=2)                #plots trace on raw image frame
    plt.title('Trace marked with a line, visit {0}, orbit {1}'.format(ancil.visnum[i], ancil.orbnum[i]))
    if ancil.trace_image_output:
        if not os.path.isdir('diagnostic_plots/trace/'):
            os.makedirs('diagnostic_plots/trace/')
        plt.savefig('diagnostic_plots/trace/{0}.png'.format(i))
    if ancil.trace_image_show:
        plt.show()
    plt.close()
    return [trace_i, trace_f]

###############################################################################################################################################################
###############################################################################################################################################################

#STEP 0: User-set parameters and constants

if ancil.output == True:
    dirname = "extracted_lc/" + datetime.strftime(datetime.now(), '%Y-%m-%d_%H:%M')
    if not os.path.exists(dirname): os.makedirs(dirname)

    table_white = QTable(names=('phase', 'spec_opt', 'var_opt', 'spec_box', 'var_box', 't_mjd', 'nvisit','norbit','scan','t_visit','t_orbit', 't_bjd'))
    table_spec = QTable(names=('t_mjd', 'phase', 'spec_opt', 'var_opt', 'template_waves', 'nvisit', 'norbit', 'scan', 't_bjd'))
    table_diagnostics = QTable(names=('nspectra', 't_mjd', 'numoutliers', 'skymedian', 'shift', "# nans"))

    copyfile("./config/obs_par.yaml", dirname+"/obs_par.yaml")        #stores obs_par.txt
    copyfile("./config/fit_par.txt", dirname+"/fit_par.txt")        #stores obs_par.txt


nspectra = 0                                                #iterator variable to track number of spectra reduced
for i, f in enumerate(ancil.files[0:2]):      #[0:104]                   #ancil.files only contains the spectra not the di
    print("\nProgress: {0}/{1}".format(i, len(ancil.files)))
    print("Filename: {0}".format(f))
    d = fits.open(f)                                        #opens the file
    scan = ancil.scans[i]


    #print(d[1].data)
    plt.imshow(d[1].data, origin = 'lower',  vmin=0, vmax=300)
    plt.colorbar()
    plt.tight_layout()
    plt.title('Background, visit {0}, orbit {1}'.format(ancil.visnum[i], ancil.orbnum[i]))
    if ancil.all_image_output:
        if not os.path.isdir('diagnostic_plots/pics/'):
            os.makedirs('diagnostic_plots/pics/')
        plt.savefig('diagnostic_plots/pics/background{0}.png'.format(i))
    if ancil.all_image_show:
        plt.show()
    plt.close()

    tmp = sum((d[1].data)[i] for i in range(len(d[1].data)))
    plt.plot(range(522), tmp)
    plt.show()
    np.savetxt('L89-59-firstexp.txt', list(zip(range(522), tmp)))

    #ancil.visnum, ancil.orbnum = ancil.visnum[i], ancil.orbnum[i]

    visnum, orbnum = ancil.visnum[i], ancil.orbnum[i]     #finds visit number and orbit number (starts at 0)
    print('visit, orbit:', (visnum, orbnum))

    if ancil.trace_image_output or ancil.trace_image_show:
        plot_trace(d, ancil)


    offset = obs_par['offset']
    cmin = int(ancil.refpix[orbnum,2] + ancil.POSTARG1/ancil.platescale) + ancil.BEAMA_i + ancil.LTV1 + offset                      #determines left column for extraction (beginning of the trace)
    cmax = min(int(ancil.refpix[orbnum,2] + ancil.POSTARG1/ancil.platescale) + ancil.BEAMA_f + ancil.LTV1 - offset, ancil.subarray_size)     #right column (end of trace, or edge of detector)
    rmin, rmax = int(obs_par['rmin']), int(obs_par['rmax'])                     #top and bottom row for extraction (specified in obs_par.txt)

    skycmin, skycmax, skyrmin, skyrmax = ancil.skycmin, ancil.skycmax, ancil.skyrmin, ancil.skyrmax

    D = np.zeros_like(d[1].data[rmin:rmax,cmin:cmax])                        #array to store the background-subtracted data
    outlier_array = np.zeros_like(d[1].data[rmin:rmax,cmin:cmax])                    #array used to determine which pixel is the biggest outlier    
    M = np.ones_like(d[1].data[rmin:rmax, cmin:cmax])                        #mask for bad pixels



    ancil.wave_grid = get_wave_grid(ancil)  # gets grid of wavelength solutions for each orbit and row
    flatfield = get_flatfield(ancil)


    bpix = d[3].data[rmin:rmax,cmin:cmax]
    badpixind =  (bpix==4)|(bpix==512)|(flatfield[orbnum][rmin:rmax, cmin:cmax] == -1.)    #selects bad pixels
    M[badpixind] = 0.0                                        #initializes bad pixel mask    

    spec_box = np.zeros(cmax - cmin)                                #box extracted standard spectrum
    spec_opt = np.zeros(cmax - cmin)                                #optimally extracted spectrum
    var_box = np.zeros(cmax - cmin)                                #box spectrum variance    
    var_opt = np.zeros(cmax - cmin)                                #optimal spectrum variance

    print((cmin, cmax))

    plt.imshow(d[1].data, vmin=0, vmax=500, origin='lower')
    plt.plot([cmin, cmin, cmax, cmax, cmin], [rmin, rmax, rmax, rmin, rmin], lw=3, c='r', alpha=0.85)
    plt.plot([skycmin, skycmin, skycmax, skycmax, skycmin], [skyrmin, skyrmax, skyrmax, skyrmin, skyrmin], lw=2, c='w',
             ls='--', alpha=0.85)
    plt.colorbar()
    plt.title('Boxes, visit {0}, orbit {1}'.format(ancil.visnum[i], ancil.orbnum[i]))
    if ancil.all_image_output:
        if not os.path.isdir('diagnostic_plots/pics/'):
            os.makedirs('diagnostic_plots/pics/')
        plt.savefig('diagnostic_plots/pics/boxes{0}.png'.format(i))
    if ancil.all_image_show:
        plt.show()
    plt.close()


    #########################################################################################################################################################
    # loops over up-the-ramp-samples (skipping first two very short exposures); gets all needed input for optextr routine                    #
    #########################################################################################################################################################

    for ii in range(d[0].header['nsamp']-2):
    #for ii in range(1):
        print("Progress (up-the-ramp-samples): {0}/{1}".format(ii, d[0].header['nsamp']-2-1))
        diff = d[ii*5 + 1].data[rmin:rmax,cmin:cmax] - d[ii*5 + 6].data[rmin:rmax,cmin:cmax]    #creates image that is the difference between successive scans

        plt.imshow(diff, vmin=0, vmax=500, origin = 'lower')
        plt.colorbar()
        plt.tight_layout()
        plt.title('UpTheRamp {0}-{1}, visit {2}, orbit {3}'.format(i,ii,ancil.visnum[i], ancil.orbnum[i]))
        if ancil.all_image_output:
            if not os.path.isdir('diagnostic_plots/pics/'):
                os.makedirs('diagnostic_plots/pics/')
            plt.savefig('diagnostic_plots/pics/uptheramp{0}-{1}.png'.format(i,ii))
        if ancil.all_image_show:
            plt.show()
        plt.close()

        diff = diff/flatfield[orbnum][rmin:rmax, cmin:cmax]                               #flatfields the differenced image

        idx = np.argmax(scipy.signal.medfilt(np.sum(diff, axis = 1),3))                         #computes spatial index of peak counts in the difference image
        #estimates sky background and variance 


        #print(d[1].data.shape)


        fullframe_diff = d[ii*5 + 1].data - d[ii*5 + 6].data                                       #fullframe difference between successive scans
        #fullframe_diff = d[1].data 

        skymedian = np.median(fullframe_diff[skyrmin:skyrmax,skycmin:skycmax])                    #estimates the background counts
        skyvar = median_abs_dev(fullframe_diff[skyrmin:skyrmax,skycmin:skycmax].flatten())            #variance for the background count estimate

        diff = diff - skymedian                                    #subtracts the background


        #interpolation to correct for bad pixels and the fact that the wavelength solution changes row by row 
        """for jj in range(rmax):
            goodidx = M[jj,:] == 1.0                            #selects good pixels
            diff[jj,:] = np.interp(ancil.wave_grid[0, int(ancil.refpix[0,1]) + ancil.LTV1, \
                cmin:cmax], ancil.wave_grid[0, jj, cmin:cmax][goodidx], diff[jj,goodidx])    #LK interpolation 8/18"""


        spectrum = diff[max(idx-ancil.window, 0):min(idx+ancil.window, rmax),:]        #selects postage stamp centered around spectrum

        #stores median of column that has spectrum on it (+/- 5 pix from center) for ii = 0

        err = np.zeros_like(spectrum) + float(obs_par['rdnoise'])**2 + skyvar
        var = abs(spectrum) + float(obs_par['rdnoise'])**2 +skyvar                #variance estimate: Poisson noise from photon counts (first term)  + readnoise (factor of 2 for differencing) + skyvar
        spec_box_0 = spectrum.sum(axis = 0)                            #initial box-extracted spectrum 
        var_box_0 = var.sum(axis = 0)                                #initial variance guess
        newM = np.ones_like(spectrum)  # not masking any pixels because we interpolated over them
        #print(M.shape)
        #print(newM.shape)
        if obs_par['opt_extract']==True: [f_opt_0, var_opt_0, numoutliers] = optextr.optextr(spectrum, err, spec_box_0, var_box_0, newM, ancil.nsmooth, ancil.sig_cut, ancil.diagnostics)
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

    phase = 0

    shift = 0.
    #corrects for wavelength drift over time
    if obs_par['correct_wave_shift'] == True:
        if nspectra == 0:
            template_waves = ancil.wave_grid[0, int(ancil.refpix[0,1]) + ancil.LTV1, cmin:cmax]             #LK interpolation 8/18 #use stellar model instead
            #template_waves -= 70.               #LK hack to get the wavelength solution right   
            print("LK adjusting wavelength solution BY HAND!!!!")
            shift = 0.
            template_spectrum = spec_opt                                #makes the first exposure the template spectrum #stellar model * grism throughput
            best_spec = spec_opt
            best_var = var_opt
        else:
            #shifts spectrum so it matches the template
            [best_spec, best_var, shift] = interpolate_spectrum(spec_opt, np.sqrt(var_opt), template_spectrum, template_waves)
            best_var = best_var**2
                
            spec_opt = best_spec                                    #saves the interpolated spectrum
            var_opt = best_var


    wavelengthsolutionoffset = 105.

    #plots wavelength solution compared to a stellar template
    """plt.plot(template_waves - wavelengthsolutionoffset, best_spec/np.max(best_spec), label = 'LK spectrum')
    np.save("temp_spectrum", [template_waves, best_spec])

    synth = np.load("syntehtic_spec.npz.npy")
    plt.plot(synth[0], 1.05*synth[1]/np.max(synth[1]), label = 'template')
    plt.legend()
    plt.show()"""


    #print(phase, sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time, visnum, orbnum, scan)

   # print(phase[0], sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time[0], visnum, orbnum, scan)
    if ancil.output == True:
        table_white.add_row([phase, sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), ancil.t_mjd[i], visnum, orbnum, scan, ancil.t_visit[i], ancil.t_orbit[i], ancil.t_bjd[i]])
        #ascii.write(table_white, dirname+'/lc_white.txt', format='ecsv', overwrite=True)
        #print(phase[0], sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time[0], visnum, orbnum, scan, t_visiti, t_orbiti, file=whitefile)
        n = len(spec_opt)
        for ii in arange(n):
            #print(time[0], phase[0], spec_opt[ii], var_opt[ii], template_waves[ii], visnum, orbnum, scan, file=specfile)
            #print(time[0], phase[0], spec_opt[ii], var_opt[ii], template_waves[ii] - wavelengthsolutionoffset, visnum, orbnum, scan, file=specfile)
            table_spec.add_row([ancil.t_mjd[i], phase, spec_opt[ii], var_opt[ii], template_waves[ii] - wavelengthsolutionoffset, visnum, orbnum, scan, ancil.t_bjd[i]])
        #print(nspectra, time[0], numoutliers, skymedian, shift, file=diagnosticsfile)
        table_diagnostics.add_row([nspectra, ancil.t_mjd[i], numoutliers, skymedian, shift, sum(np.isnan(spec_opt))])

    nspectra += 1
    #if nspectra%10 == 0: print("Extraction", '{0:1f}'.format(float(nspectra)/float(len(ancil.files))*100.), "% complete, time elapsed (min) =", '{0:0.1f}'.format(clock()/60.))


    print("sky background!!", skymedian, visnum)
    #print("length of spectrum", len(spec_opt))
    print("# nans", sum(np.isnan(spec_opt)))
    print(nspectra, ancil.t_bjd[i], sum(spec_opt), np.sum(spec_box), visnum, f, shift)

if ancil.output == True:
    ascii.write(table_white, dirname+'/lc_white.txt', format='ecsv', overwrite=True)
    ascii.write(table_spec, dirname+'/lc_spec.txt', format='ecsv', overwrite=True)
    ascii.write(table_diagnostics, dirname+'/diagnostics.txt', format='ecsv', overwrite=True)
    #specfile.close()
    #whitefile.close()
    #diagnosticsfile.close()

print("I made a change to how the wavelength interpolation is being done")
print("it's commented as LK 8/18")
print("Should check at some point to see if it's right")

