import numpy as np
from astropy.io import ascii, fits
import math
import numpy as np
from importlib import reload
import multiprocessing as mp
import os
from tqdm import tqdm
import numpy.ma as ma


def computeRMS(data, maxnbins=None, binstep=1, isrmserr=False):
    # bin data into multiple bin sizes
    npts    = data.size
    if maxnbins == None:
            maxnbins = npts/10.
    binsz   = np.arange(1, maxnbins+binstep, step=binstep)
    nbins   = np.zeros(binsz.size)
    rms     = np.zeros(binsz.size)
    rmserr  = np.zeros(binsz.size)
    for i in range(binsz.size):
            nbins[i] = int(np.floor(data.size/binsz[i]))
            bindata   = np.zeros(nbins[i], dtype=float)
# bin data
# ADDED INTEGER CONVERSION, mh 01/21/12
            for j in range(int(nbins[i])):
                    bindata[j] = data[j*binsz[i]:(j+1)*binsz[i]].mean()
            # get rms
            rms[i]    = np.sqrt(np.mean(bindata**2))
            rmserr[i] = rms[i]/np.sqrt(2.*int(nbins[i]))
    # expected for white noise (WINN 2008, PONT 2006)
    stderr = (data.std()/np.sqrt(binsz))*np.sqrt(nbins/(nbins - 1.))
    if isrmserr == True:
            return rms, stderr, binsz, rmserr
    else:
            return rms, stderr, binsz


def quantile(x, q):
    return np.percentile(x, [100. * qi for qi in q])

def weighted_mean(data, err):            #calculates the weighted mean for data points data with std devs. err
    ind = err != 0.0
    weights = 1.0/err[ind]**2
    mu = np.sum(data[ind]*weights)/np.sum(weights)
    var = 1.0/np.sum(weights)
    return [mu, np.sqrt(var)]                




def readfiles(meta):
    """
    Reads in the files saved in topdir + datadir and saves them into a list

    Args:
        meta: metadata object

    Returns:
        meta: metadata object but adds segment_list to metadata containing the sorted data fits files
    """
    meta.segment_list = []
    for fname in os.listdir(str(meta.path)):
        if fname.endswith(meta.suffix + '.fits'):
            meta.segment_list.append(str(meta.path) +'/'+ fname)
    return meta


def ancil(meta, s10=False, s20=False):

    filelist = ascii.read(meta.workdir + '/filelist.txt')

    aa = filelist['iorbit'].data
    bb = np.diff(aa)
    meta.norbit, meta.nvisit = len(np.insert(aa[1:][np.where(bb != 0)], 0, aa[0])), len(set(filelist['ivisit'].data))

    meta.mask_sp = np.array([i[0] for i in filelist['filter/grism']]) == 'G'
    meta.mask_di = np.array([i[0] for i in filelist['filter/grism']]) == 'F'
    meta.files_sp = [meta.path + '/' + i for i in filelist['filenames'][meta.mask_sp].data]
    meta.files_di = [meta.path + '/' + i for i in filelist['filenames'][meta.mask_di].data]

    f = fits.open(meta.files_sp[0])

    meta.ra = f[0].header['ra_targ'] * math.pi / 180.0  # stores right ascension
    meta.dec = f[0].header['dec_targ'] * math.pi / 180.0  # stores declination

    meta.coordtable = []  # table of spacecraft coordinates
    for i in range(max(filelist['ivisit']) + 1): meta.coordtable.append(
        meta.workdir + '/ancil/horizons/' + "/horizons_results_v" + str(i) + ".txt")

    ###
    # 03
    meta.filter = filelist['filter/grism'][meta.mask_di][0]
    meta.grism = filelist['filter/grism'][meta.mask_sp][0]

    meta.scans_sp = filelist['scan'][meta.mask_sp].data
    meta.iorbit_sp = filelist['iorbit'][meta.mask_sp].data
    meta.ivisit_sp = filelist['ivisit'][meta.mask_sp].data

    c = np.zeros(len(meta.iorbit_sp), dtype=int)
    for i in range(len(c) - 1):
        cc = np.diff(meta.iorbit_sp)
        if cc[i] == 0:
            c[i + 1] = c[i]
        else:
            c[i + 1] = c[i] + 1
    meta.iorbit_sp_com  = c

    meta.iorbit_di = filelist['iorbit'][meta.mask_di].data
    meta.ivisit_di = filelist['ivisit'][meta.mask_di].data
    meta.t_mjd_sp = filelist['t_mjd'][meta.mask_sp].data

    meta.t_orbit_sp = filelist['t_orbit'][meta.mask_sp].data
    meta.t_visit_sp = filelist['t_visit'][meta.mask_sp].data

    meta.platescale = 0.13  # IR detector has plate scale of 0.13 arcsec/pixel

    meta.POSTARG1 = f[0].header['POSTARG1']  # x-coordinate of the observer requested target offset
    meta.POSTARG2 = f[0].header['POSTARG2']  # y-coordinate of the observer requested target offset
    meta.LTV1 = int(f[1].header['LTV1'])
    meta.LTV2 = int(f[1].header['LTV2'])

    meta.subarray_size = f[1].header['SIZAXIS1']  # size of subarray

    if s10:
        if 't_bjd' in filelist.keys():
            meta.t_bjd_sp = filelist['t_bjd'][meta.mask_sp].data
        if 't_bjd' in filelist.keys():
            meta.t_bjd_di = filelist['t_bjd'][meta.mask_di].data

    if s20:
        meta.refpix = np.genfromtxt(meta.workdir + "/xrefyref.txt")  # reads in reference pixels for each visit and sorts them by time
        #idx = np.argsort(refpix[:, 0])  # sort by time
        #meta.refpix = refpix[idx]  # reference pixels from direct image

    return meta


def get_wave_grid(meta):

    if meta.grism == 'G102':
        from ..lib import geometry102 as geo
    elif meta.grism == 'G141':
        from ..lib import geometry as geo
    else:
        print('Error: GRISM in obs_par.cf is neither G102 nor G141!')

    wave_grid = np.empty((meta.norbit*meta.nvisit, meta.subarray_size, meta.subarray_size))

    #calculates wavelength solution row by row for each orbit
    for i in range(meta.norbit):
        for j in range(meta.subarray_size):
            disp_solution = geo.dispersion(meta.refpix[i,1], -meta.LTV2+j)
            delx = 0.5 + np.arange(meta.subarray_size) - (meta.refpix[i,2] + meta.LTV1 + meta.POSTARG1/meta.platescale)
            wave_grid[i, j, :] = disp_solution[0] + delx*disp_solution[1]

    return wave_grid


def get_flatfield(meta):                    #function that flatfields a data array D, which starts at [minr, cmin] of hdu[1].data
    flat = fits.open(meta.flat)                #reads in flatfield cube
    WMIN = flat[0].header['WMIN']                #constants for computing flatfield coefficients
    WMAX = flat[0].header['WMAX']

    a0 = flat[0].data[-meta.LTV1:-meta.LTV1+meta.subarray_size, -meta.LTV2:-meta.LTV2+meta.subarray_size]
    a1 = flat[1].data[-meta.LTV1:-meta.LTV1+meta.subarray_size, -meta.LTV2:-meta.LTV2+meta.subarray_size]
    a2 = flat[2].data[-meta.LTV1:-meta.LTV1+meta.subarray_size, -meta.LTV2:-meta.LTV2+meta.subarray_size]

    flatfield = []
    for i in range(meta.norbit*meta.nvisit):
        wave = meta.wave_grid[i, :]
        x = (wave - WMIN)/(WMAX-WMIN)
        flatfield.append(a0+a1*x+a2*x**2)
        flatfield[i][flatfield[i] < 0.5] = -1.        #sets flatfield pixels below 0.5 to -1 so they can be masked
    return flatfield







def median_abs_dev(vec):
    med = ma.median(vec)
    return ma.median(abs(vec - med))















def residuals(params, template_waves, template, spectrum, error):
    shift, scale = params
    fit = scale*np.interp(template_waves, template_waves-shift, spectrum)
    x = (template-fit)/error
    return (template-fit)/error

from scipy.interpolate import interp1d

def residuals2(params, x1, y1, x2, y2):
    a, b, c = params
    x1=np.array(x1)
    x2=np.array(x2)
    y1=np.array(y1)
    y2=np.array(y2)

    f = interp1d(x1, y1, kind='cubic')
    fit = f(a+b*x2)*c

    return fit - y2

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



#def getphase(t):
#    phase = (t - t0)/period
#    return phase - int(phase)









