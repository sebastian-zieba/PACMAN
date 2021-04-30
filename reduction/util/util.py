import numpy as np

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

