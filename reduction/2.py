import sys
sys.path.insert(0, './util')
import os, glob, pyfits, scipy, numpy
from numpy import *
import scipy.signal
import matplotlib.pyplot as plt
from pylab import *
import numpy.ma as ma
import geometry as geo
import suntimecorr
import optextr 
import copy
from scipy.optimize import leastsq
from astropy.io import ascii
from time import clock
from datetime import datetime
from shutil import copyfile
import pickle

class AncillaryData:
    """
    doc
    """
    def __init__(self, obs_par, fit_par, tstart, refpix):
        if not os.path.isfile("config/files.pic"):
            files = []    
            flist = glob.glob(os.path.join(obs_par['path'], "*ima.fits"))        #gets list of filenames in directory
            for f in flist: 
                d = pyfits.open(f)
                if d[0].header['filter'] == obs_par['GRISM']: files.append(f)
            pickle.dump(files, open("config/files.pic", "wb"))
            self.files = files
        else: self.files = pickle.load(open("config/files.pic", "rb"))

        f = pyfits.open(self.files[0])                

        self.flat = obs_par['flat']
        self.LTV1 = int(f[1].header['LTV1'])
        self.LTV2 = int(f[1].header['LTV2'])

        self.POSTARG1 = int(f[0].header['POSTARG1'])
        self.POSTARG2 = int(f[0].header['POSTARG2'])
        self.platescale = 0.13
        
        self.subarray_size = int(f[1].header['SIZAXIS1'])                #size of subarray
        self.BEAMA_i = int(obs_par['BEAMA_i'])                    #start of first order trace 
        self.BEAMA_f = int(obs_par['BEAMA_f'])                    #end of first order trace 
        self.npix = self.BEAMA_f - self.BEAMA_i                    #length of trace
        self.grism = obs_par['GRISM']


        self.plot_trace = convert_to_bool(obs_par['plot_trace'])     
        self.diagnostics = convert_to_bool(obs_par['diagnostics'])
        self.output = convert_to_bool(obs_par['output'])
        if self.output == False: print("NOTE: output is set to False!")

        self.sig_cut = float(obs_par['sig_cut'])
        self.nsmooth = int(obs_par['nsmooth'])
        self.window = int(obs_par['window'])                    #window outside of which the background is masked
        self.nvisit = int(obs_par['nvisit'])
        self.norb = int(obs_par['norb'])

        self.t0 = float(fit_par['t0'])
        self.period = float(fit_par['per'])

        self.one_di_per_visit = convert_to_bool(obs_par['one_di_per_visit'])    
        if self.one_di_per_visit==True: norb = nvisit
        
        self.ra = f[0].header['ra_targ']*math.pi/180.0                #stores right ascension    
        self.dec = f[0].header['dec_targ']*math.pi/180.0            #stores declination    
        
        self.expstart = f[0].header['expstart']                    #exposure start time 
        self.exptime = f[0].header['exptime']                    #exposure time [seconds]

        idx = np.argsort(refpix[:,0])
        self.refpix = refpix[idx]                        #reference pixels from direct image

        self.torbstart = self.refpix[:,0]                    #start times for each orbit
        self.torbstart = np.append(self.torbstart, 1.0e10)            #appends large value to make get_orbnum routine work

        self.tstart = tstart

        self.visnum = None
        self.orbnum = None
        
        self.wavegrid = None
    
        self.coordtable = []                            #table of spacecraft coordinates
        for i in range(self.nvisit): self.coordtable.append("bjd_conversion/horizons_results_v"+str(i)+".txt")    


def make_dict(table):
    return {x['parameter']: x['value'] for x in table}

def convert_to_bool(str):
    if str=="True": return True
    elif str=="False": return False
    else: return "String not equal to True or False"

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


def str_to_num(str):
    'Return value of numeric literal string or ValueError exception'
 
    # Handle '0'
    if str == '0': return 0
 
    # Int/Float/Complex
    try:
        return int(str)
    except ValueError:
        pass
    try:
        return float(str)
    except ValueError:
        pass
    return str

def read_dict(filename):
    #with open(filename,"r") as text:
    #    return dict(line.strip().split() for line in text)
    dict = {}
    with open(filename,"r") as text:
        for line in text:
            key, value = line.split()
            value = str_to_num(value)
            if value == 'True': value = True
            if value == 'False': value = False
            dict[key] = value
    return dict

def median_abs_dev(vec):
    med = ma.median(vec)
    return ma.median(abs(vec - med))

def getphase(t):
    phase = (t - t0)/period
    return phase - int(phase)

def get_flatfield(ancil):                    #function that flatfields a data array D, which starts at [minr, cmin] of hdu[1].data
    flat = pyfits.open(ancil.flat)                #reads in flatfield cube
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

#input: time, array of starting times (marking the begininning of each visit)
#returns: index corresponding to the visit number for time t (starts at 0)
def get_visnum(t, ancil):
    k = 0
    while(ancil.tstart[k] < t):
        k = k + 1
    return k-1

def get_orbnum(t, ancil):
    k = 0
    while(ancil.torbstart[k] < t):
        k = k + 1
    k = k - 1                                                    
    return k


def plot_trace(d, ancil):
    trace = geo.trace(ancil.refpix[:,1], ancil.refpix[:,2])                #determines trace coefficients
    orbnum = ancil.orbnum
    trace_i = ancil.refpix[orbnum,2] + ancil.POSTARG1/ancil.platescale + float(ancil.BEAMA_i + ancil.LTV1)          #start of trace
    trace_f = ancil.refpix[orbnum,2] + ancil.POSTARG1/ancil.platescale + float(ancil.BEAMA_f + ancil.LTV1)          #end of trace
    tracex = np.linspace(trace_i, trace_f,100)                    #x values over which to compute trace
    tracey = ancil.refpix[orbnum,1] + ancil.LTV2 + trace[0][orbnum] + \
    trace[1][orbnum]*(tracex - tracex[0])                    #y values of trace
    plt.imshow(d[1].data, origin = 'lower', vmin=0, vmax=40000)            #plots raw image
    plt.colorbar()
    plt.plot(tracex, tracey, color='yellow', linewidth=2)                #plots trace on raw image frame
    plt.title('Trace marked with a line')
    plt.show()
    return [trace_i, trace_f]

def get_wave_grid(ancil):
    wave_grid = np.empty((ancil.norb*ancil.nvisit, ancil.subarray_size, ancil.subarray_size))

    #calculates wavelength solution row by row for each orbit
    for i in range(ancil.norb*ancil.nvisit):
        for j in range(ancil.subarray_size):
            disp_solution = geo.dispersion(ancil.refpix[i,1], -ancil.LTV2+j)    
            delx = 0.5 + np.arange(ancil.subarray_size) - (ancil.refpix[i,2] + ancil.LTV1 + ancil.POSTARG1/ancil.platescale) 
            wave_grid[i, j, :] = disp_solution[0] + delx*disp_solution[1]     

    return wave_grid

###############################################################################################################################################################
###############################################################################################################################################################


#STEP 0: User-set parameters and constants 
obs_par = make_dict(ascii.read("config/obs_par.txt", Reader=ascii.CommentedHeader))
fit_par = make_dict(ascii.read("config/fit_par.txt", Reader=ascii.CommentedHeader))

tstart = np.genfromtxt("config/tstart.txt")                     #approximate start times for each observation (strictly before the actual start)
refpix = np.genfromtxt("config/xrefyref.txt")                    #reads in reference pixels for each visit and sorts them by time

ancil = AncillaryData(obs_par, fit_par, tstart, refpix)                #stores ancillary data

if ancil.output == True:
    dirname = "extracted_lc/" + datetime.strftime(datetime.now(), '%m_%d_%H_%M')
    if not os.path.exists(dirname): os.makedirs(dirname)

    specfile = open(dirname+'/lc_spec.txt', 'w')                #output file for spectroscopic lc
    whitefile = open(dirname+'/lc_white.txt', 'w')                #output file for white light curve
    diagnosticsfile = open(dirname+'/diagnostics.txt', 'w')

    copyfile("./config/obs_par.txt", dirname+"/obs_par.txt")        #stores obs_par.txt

ancil.wave_grid = get_wave_grid(ancil)                        #gets grid of wavelength solutions for each orbit and row
flatfield = get_flatfield(ancil)

nspectra = 0                                                #iterator variable to track number of spectra reduced
for f in ancil.files:
    d = pyfits.open(f)                                        #opens the file

    t = d[0].header['expstart']
    ancil.visnum, ancil.orbnum = get_visnum(t, ancil), get_orbnum(t, ancil)                #finds visit number and orbit number (starts at 0)

    if ancil.visnum < 0:
        print("weird, visnum <0")
        continue

    offset = int(obs_par['offset'])
    cmin = int(ancil.refpix[ancil.orbnum,2] + ancil.POSTARG1/ancil.platescale) + ancil.BEAMA_i + ancil.LTV1 + offset                      #determines left column for extraction (beginning of the trace)
    cmax = min(int(ancil.refpix[ancil.orbnum,2] + ancil.POSTARG1/ancil.platescale) + ancil.BEAMA_f + ancil.LTV1 - offset, ancil.subarray_size)     #right column (end of trace, or edge of detector)
    
    rmin, rmax = int(obs_par['rmin']), int(obs_par['rmax'])                     #top and bottom row for extraction (specified in obs_par.txt)

    D = np.zeros_like(d[1].data[rmin:rmax,cmin:cmax])                        #array to store the background-subtracted data
    D_nointerp = np.zeros_like(d[1].data[rmin:rmax,cmin:cmax])                    #array to store the background-subtracted data w/o interpolating over different wavelengths in every row
    outlier_array = np.zeros_like(d[1].data[rmin:rmax,cmin:cmax])                    #array used to determine which pixel is the biggest outlier    
    M = np.ones_like(d[1].data[rmin:rmax, cmin:cmax])                        #mask for bad pixels


    bpix = d[3].data[rmin:rmax,cmin:cmax]
    badpixind =  (bpix==4)|(bpix==512)|(flatfield[ancil.orbnum][rmin:rmax, cmin:cmax] == -1.)    #selects bad pixels 
    M[badpixind] = 0.0                                        #initializes bad pixel mask    

    spec_box = np.zeros(cmax - cmin)                                #box extracted standard spectrum
    spec_opt = np.zeros(cmax - cmin)                                #optimally extracted spectrum
    var_box = np.zeros(cmax - cmin)                                #box spectrum variance    
    var_opt = np.zeros(cmax - cmin)                                #optimal spectrum variance

    scan = 0                                            #sets scan direction
    if d[0].header['postarg2'] < 0: scan = 1

    #########################################################################################################################################################
    # loops over up-the-ramp-samples (skipping first two very short exposures); gets all needed input for optextr routine                    #
    #########################################################################################################################################################

    for ii in range(d[0].header['nsamp']-2):    
    #for ii in range(1):
        #print(ii)
        diff = d[ii*5 + 1].data[rmin:rmax,cmin:cmax] - d[ii*5 + 6].data[rmin:rmax,cmin:cmax]    #creates image that is the difference between successive scans
        
        #diff = d[1].data[rmin:rmax,cmin:cmax] 

        diff = diff/flatfield[ancil.orbnum][rmin:rmax, cmin:cmax]                               #flatfields the differenced image

        idx = np.argmax(scipy.signal.medfilt(np.sum(diff, axis = 1),3))                         #computes spatial index of peak counts in the difference image

        #estimates sky background and variance 
        [skyrmin, skyrmax] = [10, 70]
        [skycmin, skycmax] = [400, 500]
        
        fullframe_diff = d[ii*5 + 1].data - d[ii*5 + 6].data                                       #fullframe difference between successive scans
        #fullframe_diff = d[1].data 

        skymedian = np.median(fullframe_diff[skyrmin:skyrmax,skycmin:skycmax])                    #estimates the background counts
        skyvar = median_abs_dev(fullframe_diff[skyrmin:skyrmax,skycmin:skycmax].flatten())            #variance for the background count estimate

        diff = diff - skymedian                                    #subtracts the background

        D_nointerp = D_nointerp+diff                                #sums up the diffs of the spectrum (with no interpolation for wavelength solution changing)


        #interpolation to correct for bad pixels and the fact that the wavelength solution changes row by row 
        """for jj in range(rmax):
            goodidx = M[jj,:] == 1.0                            #selects good pixels
            diff[jj,:] = np.interp(ancil.wave_grid[0, int(ancil.refpix[0,1]) + ancil.LTV1, \
                cmin:cmax], ancil.wave_grid[0, jj, cmin:cmax][goodidx], diff[jj,goodidx])    #LK interpolation 8/18"""


        spectrum = diff[max(idx-ancil.window, 0):min(idx+ancil.window, rmax),:]        #selects postage stamp centered around spectrum 

        #stores median of column that has spectrum on it (+/- 5 pix from center) for ii = 0
        if ii == 0:
            cRate = np.median(diff[max(idx-5, 0):min(idx+5, rmax),:], axis = 0)        #selects postage stamp centered around spectrum 
            cRate /= ancil.exptime                                                     #calculates average counts over entire exposure

        err = np.zeros_like(spectrum) + float(obs_par['rdnoise'])**2 + skyvar
        var = abs(spectrum) + float(obs_par['rdnoise'])**2 +skyvar                #variance estimate: Poisson noise from photon counts (first term)  + readnoise (factor of 2 for differencing) + skyvar
        spec_box_0 = spectrum.sum(axis = 0)                            #initial box-extracted spectrum 
        var_box_0 = var.sum(axis = 0)                                #initial variance guess
    
        newM = np.ones_like(spectrum)                                #not masking any pixels because we interpolated over them 

        if convert_to_bool(obs_par['opt_extract'])==True: [f_opt_0, var_opt_0, numoutliers] = optextr.optextr(spectrum, err, spec_box_0, var_box_0, newM, ancil.nsmooth, ancil.sig_cut, ancil.diagnostics)                            
        else: [f_opt, var_opt] = [spec_box_0,var_box_0]

        #sums up spectra and variance for all the differenced images
        spec_opt += f_opt_0                            
        var_opt += var_opt_0
        spec_box += spec_box_0
        var_box += var_box_0
        
    ######################################################################################################################################    

    time = (d[0].header['expstart'] + d[0].header['expend'])/(2.0) + 2400000.5                    #converts time to BJD_TDB; see Eastman et al. 2010 equation 4
    time = time + (32.184)/(24.0*60.0*60.0)                                    
    print("visnum, suntimecorr", ancil.visnum, (suntimecorr.suntimecorr(ancil.ra, ancil.dec, array([time]), ancil.coordtable[ancil.visnum], verbose=False))/(60.0))
    time = time + (suntimecorr.suntimecorr(ancil.ra, ancil.dec, array([time]), ancil.coordtable[ancil.visnum], verbose=False))/(60.0*60.0*24.0)    

    phase = (time-ancil.t0)/ancil.period - math.floor((time-ancil.t0)/ancil.period)
    if phase > 0.5: phase = phase - 1.0                                    #calculates orbital phase

    shift = 0.
    #corrects for wavelength drift over time
    if convert_to_bool(obs_par['correct_wave_shift']) == True:
        if nspectra == 0:
            template_waves = ancil.wave_grid[0, int(ancil.refpix[0,1]) + ancil.LTV1, cmin:cmax]             #LK interpolation 8/18
            #template_waves -= 70.               #LK hack to get the wavelength solution right   
            print("LK adjusting wavelength solution BY HAND!!!!")
            shift = 0.
            template_spectrum = spec_opt                                #makes the first exposure the template spectrum
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


    #print(phase, sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time, ancil.visnum, ancil.orbnum, scan)
    n = len(spec_opt)    
   # print(phase[0], sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time[0], ancil.visnum, ancil.orbnum, scan)
    if ancil.output == True:
        print(phase[0], sum(spec_opt), sum(var_opt),  sum(spec_box), sum(var_box), time[0], ancil.visnum, ancil.orbnum, scan, np.mean(cRate), file=whitefile)
        for ii in arange(n):
            #print(time[0], phase[0], spec_opt[ii], var_opt[ii], template_waves[ii], ancil.visnum, ancil.orbnum, scan, cRate[ii], file=specfile)
            print(time[0], phase[0], spec_opt[ii], var_opt[ii], template_waves[ii] - wavelengthsolutionoffset, ancil.visnum, ancil.orbnum, scan, cRate[ii], file=specfile)
        print(nspectra, time[0], numoutliers, skymedian, shift, file=diagnosticsfile)

    nspectra += 1
    if nspectra%10 == 0: print("Extraction", '{0:1f}'.format(float(nspectra)/float(len(ancil.files))*100.), "% complete, time elapsed (min) =", '{0:0.1f}'.format(clock()/60.))


    print("sky background!!", skymedian, ancil.visnum)
    #print("length of spectrum", len(spec_opt))
    print("# nans", sum(np.isnan(spec_opt)))
    print(nspectra, time[0], sum(spec_opt), np.sum(spec_box), ancil.visnum, f, shift)

if ancil.output == True:
    specfile.close()
    whitefile.close()
    diagnosticsfile.close()

print("I made a change to how the wavelength interpolation is being done")
print("it's commented as LK 8/18")
print("Should check at some point to see if it's right")

