import numpy as np
import pysynphot as S
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d
import yaml

obs_par_path = "config/obs_par.yaml"
with open(obs_par_path, 'r') as file:
    obs_par = yaml.safe_load(file)


import sys
sys.path.insert(0, './util')
import ancil
import yaml
from astropy.io import ascii


ancil = ancil.AncillaryData(obs_par)

filelist_path = './config/filelist.txt'
data = ascii.read(filelist_path)


Teff, logg, MH = obs_par['Teff'], obs_par['logg'], obs_par['MH']

#def Planck(b_wvl, teff):
#    h =6.62607004e-34          #m^2/kg/s
#    c =299792458.0             #m/s
#    kb =1.38064852e-23         #m^2 kg /s^2 K
#    return (2.0*h*(c**2)/(b_wvl**5))*(1.0/( np.exp( (h*c)/(b_wvl*kb*teff) )- 1.0))

if obs_par['sm'] == 'phoenix':
    sm = 'phoenix'
elif obs_par['sm'] == 'kurutz':
    sm = 'ck04models'

if ancil.grism == 'G141':
    grism = 'g141'
elif ancil.grism == 'G102':
    grism = 'g102'


stellar_spectrum = S.Icat(sm, Teff, MH, logg)
wvl = stellar_spectrum.wave/1e4 #microns
flux = stellar_spectrum.flux*1e7/(np.pi)
flux = flux/max(flux)

x = wvl
y = flux
f = interp1d(x, y, kind='cubic')



nvisit = data['nvisit']


for vi in np.unique(nvisit):
    throughput = np.loadtxt('./ancil/throughputs_and_spectra/bandpass_v{0}.txt'.format(vi)).T

    wvl_g = throughput[0]/1e4 #microns
    flux_g = throughput[1]

    plt.plot(x, y, label='stellar spectrum')
    plt.plot(wvl_g, flux_g, label='bandpass')
    plt.plot(wvl_g, f(wvl_g)*flux_g, label='stellar spectrum * bandpass')

    plt.xlim(0.7, 2)

    #plt.xscale('log')
    plt.legend()

    plt.savefig('./ancil/throughputs_and_spectra/bandpass_times_spectrum_v{0}.png'.format(vi))
    plt.close()
    np.savetxt('./ancil/throughputs_and_spectra/bandpass_times_spectrum_v{0}.txt'.format(vi), list(zip(wvl_g, f(wvl_g)*flux_g)))

