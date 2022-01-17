import numpy as np
import matplotlib.pyplot as plt
import pysynphot as S 
# https://pysynphot.readthedocs.io/en/latest/bandpass.html#observation-mode
# https://pysynphot.readthedocs.io/en/latest/appendixb.html#wfc3
# https://pysynphot.readthedocs.io/en/latest/appendixb.html#mjd
import os
from astropy.io import ascii
from tqdm import tqdm
from astropy.time import Time


times = []

for i in np.arange(0,10,1):
    time = '201{}-01-01T00:00:00'.format(i)
    times.append(time)

for i in np.arange(0,2,1):
    time = '202{}-01-01T00:00:00'.format(i)
    times.append(time)

t = Time(times, format='isot', scale='utc')
print(t.mjd)

#grism = 'g141'
grism = 'g102'

print(t.mjd)

bp = S.ObsBandpass('wfc3,ir,{0},mjd#{1}'.format(grism, t.mjd[0]))

wvl, bp_val = bp.binset, bp(bp.binset)
np.savetxt('./bandpass_{0}.txt'.format(grism), list(zip(wvl, bp_val)), header='wavelength (angstrom) \t throughput')

