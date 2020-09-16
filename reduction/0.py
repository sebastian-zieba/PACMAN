# This step creates a table of all available _ima.fits files
# The contents of this filelist.txt are:
# Filename / Filter or Grism / nvisit / norbit / Time since first exposure in orbit / Time since first exposure in visit

import sys

sys.path.insert(0, './util')

from pylab import *
import gaussfitter, ancil, os, glob
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
prevrun = False  # os.path.exists('config/filelist.txt')
print(prevrun)

if not prevrun:  # if no previous run was done
    expstart = np.zeros(len(files))
    filtergrism = np.zeros(len(files), dtype=object)

    # Will create a table with the properties of all _ima.fits files at the first run
    for i in enumerate(files):
        ima = fits.open(i[1])
        filtergrism[i[0]] = str(ima[0].header['filter'])
        expstart[i[0]] = ima[0].header['expstart']
        ima.close()

    torder = np.argsort(expstart)
    files = np.array([i.split('/')[-1] for i in files])[torder]  # files are now chronologically sorted
    times = expstart[torder]
    print(times)

    # Identify orbits and visits
    norbits = []
    nvisits = []
    norbit = 0
    nvisit = 0
    tos = []
    tvs = []
    first_exp_orbit = 0
    first_exp_visit = 0

    times_diff = np.insert(np.diff(times), 0, 0)

    for i in enumerate(times):
        # if two exposures arent in the same orbit and more than an orbital period apart -> not subsequent orbits but a new visit
        if times_diff[i[0]] * 24 * 60 > 100:
            first_exp_visit = i[0]
            first_exp_orbit = i[0]
            norbit = 0
            nvisit += 1
        # if two exposure more than 10 min apart but less than an orbital period -> subsequent orbits
        elif 10 < times_diff[i[0]] * 24 * 60 <= 100:
            first_exp_orbit = i[0]
            norbit += 1
        # else: two exposures less than 10mins apart -> same orbit
        norbits.append(norbit)
        nvisits.append(nvisit)
        to = (i[1] - times[first_exp_orbit]) * 24 * 60  # time since first exposure in orbit
        tv = (i[1] - times[first_exp_visit]) * 24 * 60  # time since first exposure in orbit
        tos.append(to)
        tvs.append(tv)

    t = QTable([files, filtergrism[torder], nvisits, norbits, tos, tvs],
               names=('filenames', 'filter/grism', 'nvisit', 'norbit', 't_orbit', 't_visit'))

    ascii.write(t, 'config/filelist.txt', format='ecsv')
