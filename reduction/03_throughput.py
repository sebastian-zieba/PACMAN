import numpy as np
import matplotlib.pyplot as plt
import pysynphot as S 

#https://pysynphot.readthedocs.io/en/latest/bandpass.html#observation-mode
#https://pysynphot.readthedocs.io/en/latest/appendixb.html#wfc3
#https://pysynphot.readthedocs.io/en/latest/appendixb.html#mjd

import sys
sys.path.insert(0, './util')
import ancil
import yaml
from astropy.io import ascii


#import obs_par and fit_par
obs_par_path = "config/obs_par.yaml"
with open(obs_par_path, 'r') as file:
    obs_par = yaml.safe_load(file)
ancil = ancil.AncillaryData(obs_par)

filelist_path = './config/filelist.txt'
data = ascii.read(filelist_path)

nvisit = data['nvisit']
t_mjd = data['t_mjd']

#store time when every visit started
t_mjd_visit_starts = t_mjd[np.unique(nvisit, return_index=True)[1]]


if ancil.grism == 'G141':
    grism = 'g141'
elif ancil.grism == 'G102':
    grism = 'g102'


for i, t_mjd_visit_start in enumerate(t_mjd_visit_starts):
    bp = S.ObsBandpass('wfc3,ir,{0},mjd#{1}'.format(grism, t_mjd_visit_start))
    wvl, bp_val = bp.binset, bp(bp.binset)

    np.savetxt('./ancil/throughputs_and_spectra/bandpass_v{0}.txt'.format(i), list(zip(wvl, bp_val)))

    plt.plot(wvl, bp_val, c='C0')
    plt.xlabel('Angstrom')
    plt.ylabel('Throughput')
    plt.title('{0}_v{1}'.format(grism, i))
    plt.savefig('./ancil/throughputs_and_spectra/bandpass_v{0}.png'.format(i))
    plt.close()


