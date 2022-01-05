import numpy as np
import pysynphot as S
import os
from astropy.io import ascii
from scipy.interpolate import interp1d
from tqdm import tqdm
import matplotlib.pyplot as plt


wvls = np.linspace(0.2, 5, 100)/1e6

Teff   = 5000
logg   = 4.5
MH     = 0

h = 6.62607004e-34  # m^2 kg/s
c = 299792458.0  # m/s
kb = 1.38064852e-23  # m^2 kg /s^2 K


sm = 'phoenix'
# Generate Stellar spectrum
sp = S.Icat(sm, Teff, MH, logg)
print(sp.waveunits.name)
print(sp.fluxunits.name)
sp_wvl = sp.wave
sp_flux_flam = sp.flux

plt.plot(sp_wvl*1e-10, sp_flux_flam*1e10)
plt.xscale('log')
plt.xlim(100e-9, 10e-6)
plt.show()

sp.convert('Fnu')
print(sp.waveunits.name)
print(sp.fluxunits.name)
sp_wvl = sp.wave
sp_flux_fnu = sp.flux

plt.plot(sp_wvl*1e-10, sp_flux_fnu*3e8/((sp_wvl*1e-10)**2))
plt.xscale('log')
plt.xlim(100e-9, 10e-6)
plt.show()