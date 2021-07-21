# This step creates a table of all available _ima.fits files
# The contents of this filelist.txt are:
# Filename / Filter or Grism / nvisit / norbit / Time since first exposure in orbit / Time since first exposure in visit

import sys
import numpy as np
import os, glob, time
from astropy.table import QTable
from astropy.io import ascii, fits
import shutil
from ..lib import logedit
from ..lib import readECF as rd

from ..lib import util

class MetaClass:
    def __init__(self):
        # initialize Univ
        # Univ.__init__(self)
        # self.initpars(ecf)
        # self.foo = 2
        return



def run00(eventlabel):

    # Initialize metadata object
    meta = MetaClass()
    meta.eventlabel = eventlabel

    # Create directories for Stage 3 processing
    #t0 = time.time()
    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
    meta.workdir = 'run_' + datetime + '_' + meta.eventlabel
    if not os.path.exists(meta.workdir):
        os.makedirs(meta.workdir)
    if not os.path.exists(meta.workdir + "/figs"):
        os.makedirs(meta.workdir + "/figs")

    # Load Eureka! control file and store values in Event object
    ecffile = 'obs_par.ecf'
    ecf = rd.read_ecf(ecffile)
    rd.store_ecf(meta, ecf)

    # Copy ecf
    shutil.copy(ecffile, meta.workdir)
    shutil.copy('fit_par.txt', meta.workdir)

    # Create list of file segments
    meta = util.readfiles(meta)
    num_data_files = len(meta.segment_list)
    print(f'\nFound {num_data_files} data file(s) ending in {meta.suffix}.fits')
    files = meta.segment_list  # gets list of filenames in directory

    # A file called "filelist.txt" should include all _ima.fits file names and if they are a spec or di
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

    # files are chronologically sorted
    tsort = np.argsort(times)
    files = np.array([i.split('/')[-1] for i in files])[tsort]
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

    for i, timei in enumerate(times):
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
        tos[i] = (timei - times[orbit_begin_i]) * 24 * 60  # time since first exposure in orbit
        tvs[i] = (timei - times[visit_begin_i]) * 24 * 60  # time since first exposure in orbit

    table = QTable([files, filter, nvisits, norbits, times, tvs, tos, scans, exp],
               names=('filenames', 'filter/grism', 'nvisit', 'norbit', 't_mjd', 't_visit', 't_orbit',
                      'scan', 'exp'))# scan: (0: forward - lower flux, 1: reverse - higher flux, -1: postarg2=0)
    ascii.write(table, meta.workdir + '/filelist.txt', format='ecsv', overwrite=True)

    #table = QTable([], names=())
    #ascii.write(table, 'config/filelist.txt', format='ecsv', overwrite=True)
    #table.add_columns([files, filter, nvisits, norbits, times, tos, tvs, scans, exp],
    #                  names=['filenames', 'filter/grism', 'nvisit', 'norbit', 't_mjd', 't_orbit', 't_visit','scan', 'exp'])# scan: (0: forward - lower flux, 1: reverse - higher flux, -1: postarg2=0)
    #ascii.write(table, 'config/filelist.txt', format='ecsv', overwrite=True)

    print('Finished 00.py')

    return meta
