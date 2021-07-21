import numpy as np
from astropy.io import ascii, fits
import math


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





import numpy as np
from importlib import reload
import multiprocessing as mp
import os
from tqdm import tqdm

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


def ancil(meta, s20=False):

    filetable = ascii.read(meta.workdir + '/filelist.txt')
    sp_mask = np.array([i[0] for i in filetable['filter/grism']]) == 'G'
    di_mask = np.array([i[0] for i in filetable['filter/grism']]) == 'F'
    meta.files_sp = [meta.path + '/' + i for i in filetable['filenames'][sp_mask].data]



    f = fits.open(meta.files_sp[0])

    meta.ra = f[0].header['ra_targ'] * math.pi / 180.0  # stores right ascension
    meta.dec = f[0].header['dec_targ'] * math.pi / 180.0  # stores declination

    meta.coordtable = []  # table of spacecraft coordinates
    for i in range(max(filetable['nvisit']) + 1): meta.coordtable.append(
        meta.workdir + '/ancil/horizons/' + "/horizons_results_v" + str(i) + ".txt")

    ###
    # 03
    meta.filter = filetable['filter/grism'][di_mask][0]
    meta.grism = filetable['filter/grism'][sp_mask][0]

    meta.scans = filetable['scan'][sp_mask].data
    meta.orbnum = filetable['norbit'][sp_mask].data
    meta.visnum = filetable['nvisit'][sp_mask].data

    meta.t_mjd = filetable['t_mjd'][sp_mask].data

    meta.t_orbit = filetable['t_orbit'][sp_mask].data
    meta.t_visit = filetable['t_visit'][sp_mask].data

    meta.platescale = 0.13  # IR detector has plate scale of 0.13 arcsec/pixel

    meta.POSTARG1 = f[0].header['POSTARG1']  # x-coordinate of the observer requested target offset
    meta.POSTARG2 = f[0].header['POSTARG2']  # y-coordinate of the observer requested target offset
    meta.LTV1 = int(f[1].header['LTV1'])
    meta.LTV2 = int(f[1].header['LTV2'])

    meta.subarray_size = f[1].header['SIZAXIS1']  # size of subarray

    if s20 == True:
        if 't_bjd' in filetable.keys():
            meta.t_bjd = filetable['t_bjd'][sp_mask].data

        refpix = np.genfromtxt(meta.workdir + "/xrefyref.txt")  # reads in reference pixels for each visit and sorts them by time
        idx = np.argsort(refpix[:, 0])  # sort by time
        meta.refpix = refpix[idx]  # reference pixels from direct image

    return meta