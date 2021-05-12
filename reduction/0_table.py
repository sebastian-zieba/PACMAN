# This step creates a table of all available _ima.fits files
# The contents of this filelist.txt are:
# Filename / Filter or Grism / nvisit / norbit / Time since first exposure in orbit / Time since first exposure in visit

import sys
sys.path.insert(0, './util')
import ancil

import numpy as np
import os, glob
import yaml
from astropy.table import QTable
from astropy.io import ascii, fits

obs_par_path = "config/obs_par.yaml"
with open(obs_par_path, 'r') as file:
    obs_par = yaml.safe_load(file)

ancil = ancil.AncillaryData(obs_par)

files = glob.glob(os.path.join(ancil.path, "*_ima.fits"))  # gets list of filenames in directory

# If file exists, a previous runs has already been executed
# Then, a file called "filelist.txt" should include all _ima.fits file names and if they are a spec or di
#prevrun = False #os.path.exists('config/filelist.txt')
#print(prevrun)

#if not prevrun:  # if no previous run was done
times = np.zeros(len(files))
exp = np.zeros(len(files))
filter = np.zeros(len(files), dtype=object)
scans = np.zeros(len(files), dtype=int) #stores scan directions

# Will create a table with the properties of all _ima.fits files at the first run
for i, file in enumerate(files):
    ima = fits.open(file)
    #the header "filter" tells us if the observation used a Filter (-> Direct Image) or a Grism (-> Spectrum)
    filter[i] = str(ima[0].header['filter'])
    exp[i] = ima[0].header['exptime']
    times[i] = (ima[0].header['expstart'] + ima[0].header['expend'])/(2.0)#ima[0].header['expstart']
    #scan direction
    scans[i] = 0  # sets scan direction
    if ima[0].header['postarg2'] < 0: scans[i] = 1
    elif ima[0].header['postarg2'] == 0: scans[i]= -1
    ima.close()

tsort = np.argsort(times)
files = np.array([i.split('/')[-1] for i in files])[tsort]  # files are now chronologically sorted
times = times[tsort]
exp = exp[tsort]
filter = filter[tsort]
scans = scans[tsort]

# Identify orbits and visits
norbits = np.zeros(len(times), dtype=int)
nvisits = np.zeros(len(times), dtype=int)
norbit = 0
nvisit = 0
tos = np.zeros(len(times)) #time since begin of orbit
tvs = np.zeros(len(times)) #time since begin of visit
orbit_begin_i = 0 #index of first exposure in orbit
visit_begin_i = 0 #index of first exposure in visit

times_diff = np.insert(np.diff(times), 0, 0)

for i, time in enumerate(times):
    # if two exposures arent in the same orbit and more than an orbital period apart -> not subsequent orbits but a new visit
    if times_diff[i] * 24 * 60 > 100:
        orbit_begin_i = i
        visit_begin_i = i
        norbit = 0
        nvisit += 1
    # if two exposure more than 10 min apart but less than an orbital period -> subsequent orbits
    elif 10 < times_diff[i] * 24 * 60 <= 100:
        orbit_begin_i = i
        norbit += 1
    # else: two exposures less than 10mins apart -> same orbit and same visit
    norbits[i] = norbit
    nvisits[i] = nvisit
    tos[i] = (time - times[orbit_begin_i]) * 24 * 60  # time since first exposure in orbit
    tvs[i] = (time - times[visit_begin_i]) * 24 * 60  # time since first exposure in orbit

table = QTable([files, filter, nvisits, norbits, times, tvs, tos, scans, exp],
           names=('filenames', 'filter/grism', 'nvisit', 'norbit', 't_mjd', 't_visit', 't_orbit',
                  'scan', 'exp'))# scan: (0: forward - lower flux, 1: reverse - higher flux, -1: postarg2=0)
ascii.write(table, 'config/filelist.txt', format='ecsv', overwrite=True)

#table = QTable([], names=())
#ascii.write(table, 'config/filelist.txt', format='ecsv', overwrite=True)
#table.add_columns([files, filter, nvisits, norbits, times, tos, tvs, scans, exp],
#                  names=['filenames', 'filter/grism', 'nvisit', 'norbit', 't_mjd', 't_orbit', 't_visit','scan', 'exp'])# scan: (0: forward - lower flux, 1: reverse - higher flux, -1: postarg2=0)
#ascii.write(table, 'config/filelist.txt', format='ecsv', overwrite=True)


print('Finished 0.py')