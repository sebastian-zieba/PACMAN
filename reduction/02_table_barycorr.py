# Based on a perl script found on https://renenyffenegger.ch/notes/Wissenschaft/Astronomie/Ephemeriden/JPL-Horizons
# Retrieves vector data of Hubble from JPL's HORIZONS system on https://ssd.jpl.nasa.gov/horizons_batch.cgi (see Web interface on https://ssd.jpl.nasa.gov/horizons.cgi)
# Also helpful: https://github.com/kevin218/POET/blob/master/code/doc/spitzer_Horizons_README.txt

import sys
sys.path.insert(0, './util')
import os
import numpy as np
from astropy.io import ascii
#from astropy.time import Time
import urllib.request 
import suntimecorr
import ancil
import yaml
import math
from astropy.table import QTable, Column

#import obs_par and fit_par
obs_par_path = "config/obs_par.yaml"
with open(obs_par_path, 'r') as file:
	obs_par = yaml.safe_load(file)

ancil = ancil.AncillaryData(obs_par)

filelist_path = './config/filelist.txt'
data = ascii.read(filelist_path)

nvisit = data['nvisit']
t_mjd = data['t_mjd']
t_bjd = np.zeros(len(t_mjd))

for i in range(max(nvisit)+1):
	ivisit = nvisit == i
	t_jd = t_mjd[ivisit] + 2400000.5  # converts time to BJD_TDB; see Eastman et al. 2010 equation 4
	t_jd = t_jd + (32.184) / (24.0 * 60.0 * 60.0)
	t_bjd[ivisit] = t_jd + (suntimecorr.suntimecorr(ancil.ra, ancil.dec, t_jd, ancil.coordtable[i], verbose=False)) / (60.0 * 60.0 * 24.0)
	#phase = (time - ancil.t0) / ancil.period - math.floor((time - ancil.t0) / ancil.period)
	#if phase > 0.5: phase = phase - 1.0

if not any(np.array(data.keys())=='t_bjd'):
	data.add_column(Column(data = t_bjd, name='t_bjd'))
	ascii.write(data, filelist_path, format='ecsv', overwrite=True)
else:
	data.replace_column(name='t_bjd', col=Column(data = t_bjd, name='t_bjd'))
	ascii.write(data, filelist_path, format='ecsv', overwrite=True)