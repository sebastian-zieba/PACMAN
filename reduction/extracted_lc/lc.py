import numpy as np
import matplotlib.pyplot as plt
import batman
from astropy.io import ascii
from scipy.optimize import curve_fit



datafile = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/reduction/extracted_lc/2021-04-21_19:42/lc_white.txt'

data = np.loadtxt(datafile).T

t = data[5]-2400000.5
f = data[1]/np.median(data[1])


filetable = ascii.read('/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/reduction/config/filelist.txt')

#only use G141 data
#only use first visit
mask = (filetable['filter/grism'] == 'G141') & (filetable['nvisit'] == 0)


orbnum = filetable['norbit'][mask].data
visnum = filetable['nvisit'][mask].data

scan = (filetable['scan'][mask].data == 0)[26:]


torb0 = t[np.where(orbnum == 0)[0][0]]
torb1 = t[np.where(orbnum == 1)[0][0]]
torb2 = t[np.where(orbnum == 2)[0][0]]
torb3 = t[np.where(orbnum == 3)[0][0]]
torbs = np.array([torb0, torb1, torb2, torb3])
torbs = np.array([torb1, torb2, torb3])

tvis0 = t[np.where(visnum == 0)[0][0]]
tvis0 = torb1

#remove first orbit
q = t > 58888.3

t = t[q]
f = f[q]
plt.scatter(t,f, label='raw', c='r', s=10)


# model

def exponential_slope(t, a, b, c, v, si):
    to=t - np.array([torbs[[ti >= np.array(torbs) for ti in t][ii]][-1] for ii in range(len(t))])
    tv=t - tvis0
    s = np.ones(len(t))
    s[scan] = si

    params = batman.TransitParams()
    params.t0 = 58365.6708  # +2400000.5                        #time of inferior conjunction
    params.per = 2.2531  # orbital period
    params.rp = 0.0234  # planet radius (in units of stellar radii)
    params.a = 15.829  # semi-major axis (in units of stellar radii)
    params.inc = 90  # orbital inclination (in degrees)
    params.ecc = 0.  # eccentricity
    params.w = 90.  # longitude of periastron (in degrees)
    params.u = [0.1, 0.3]  # limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"  # limb darkening model
    m = batman.TransitModel(params, t)  # initializes model
    fluxb = m.light_curve(params)  # calculates light curve

    return fluxb*s*((1-np.exp(-a*to-b))*(c+v*tv))




#fit model to raw data

init_vals = [100, 1, 1, 1, 1.01]

popt, pcov = curve_fit(exponential_slope, t, f, p0=init_vals)
print('best_vals: {}'.format(popt))

plt.plot(t, exponential_slope(t, *popt), 'r--', alpha=0.4, label='model: transit*s*((1-np.exp(-a*to-b))*(c+v*tv))')
plt.scatter(t, f-exponential_slope(t, *popt)+1, label='residuals', s=20, c='k', zorder=-10)


# plot transit

tmodel = np.linspace(t[0],t[-1], 1000)
params = batman.TransitParams()
params.t0 = 1366.1701+2457000-2400000.5                        #time of inferior conjunction
params.per = 2.25314                      #orbital period
params.rp = 0.0234                     #planet radius (in units of stellar radii)
params.a = 16.2                       #semi-major axis (in units of stellar radii)
params.inc = 88.7                     #orbital inclination (in degrees)
params.ecc = 0.                      #eccentricity
params.w = 90.                       #longitude of periastron (in degrees)
params.u = [0.1, 0.3]                #limb darkening coefficients [u1, u2]
params.limb_dark = "quadratic"       #limb darkening model
m = batman.TransitModel(params, tmodel)    #initializes model
fluxb = m.light_curve(params)          #calculates light curve
plt.plot(tmodel,fluxb, label='transit (lit.)')




plt.xlim(58888.3, 58888.5)
plt.ylim(0.9975, 1.0015)

plt.legend()
plt.savefig('L98-59b_transit_visit0.jpeg', dpi=250)
plt.show()









