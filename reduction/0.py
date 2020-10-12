# This step creates a table of all available _ima.fits files
# The contents of this filelist.txt are:
# Filename / Filter or Grism / nvisit / norbit / Time since first exposure in orbit / Time since first exposure in visit

import sys
sys.path.insert(0, './util')

#from pylab import *
import numpy as np
import ancil, os, glob
from astropy.io import fits
import yaml
from astropy.table import QTable
from astropy.io import ascii

obs_par_path = "config/obs_par.yaml"
with open(obs_par_path, 'r') as file:
    obs_par = yaml.safe_load(file)

ancil = ancil.AncillaryData(obs_par)

files = glob.glob(os.path.join(ancil.path, "*_ima.fits"))  # gets list of filenames in directory

# If file exists, a previous runs has already been executed
# Then, a file called "filelist.txt" should include all _ima.fits file names and if they are a spec or di
prevrun = False #os.path.exists('config/filelist.txt')
print(prevrun)

if not prevrun:  # if no previous run was done
    expstart = np.zeros(len(files))
    filter = np.zeros(len(files), dtype=object)

    # Will create a table with the properties of all _ima.fits files at the first run
    for i, file in enumerate(files):
        ima = fits.open(file)
        #the header "filter" tells us if the observation used a Filter (-> Direct Image) or a Grism (-> Spectrum)
        filter[i] = str(ima[0].header['filter'])
        expstart[i] = ima[0].header['expstart']
        ima.close()

    tsort = np.argsort(expstart)
    files = np.array([i.split('/')[-1] for i in files])[tsort]  # files are now chronologically sorted
    times = expstart[tsort]
    filter = filter[tsort]

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

    t = QTable([files, filter, nvisits, norbits, tos, tvs],
               names=('filenames', 'filter/grism', 'nvisit', 'norbit', 't_orbit', 't_visit'))

    ascii.write(t, 'config/filelist.txt', format='ecsv')
