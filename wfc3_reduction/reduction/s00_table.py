# This step creates a table of all available _ima.fits files
# The contents of this filelist.txt are:
# Filename / Filter or Grism / nvisit / norbit / Time since first exposure in orbit / Time since first exposure in visit

import numpy as np
import os, time
import shutil
from astropy.table import QTable
from astropy.io import ascii, fits
from ..lib import readECF as rd
from ..lib import util
from ..lib import manageevent as me
from scipy.stats import rankdata
from tqdm import tqdm
# TODO Change name of readECF


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
    for i, file in enumerate(tqdm(files, desc='Reading in files and their header information')):
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
    iorbits = np.zeros(len(times), dtype=int)
    ivisits = np.zeros(len(times), dtype=int)
    iorbit = 0
    ivisit = 0
    tos = np.zeros(len(times)) #time since begin of orbit
    tvs = np.zeros(len(times)) #time since begin of visit
    iorbit_begin = 0 #index of first exposure in orbit
    ivisit_begin = 0 #index of first exposure in visit

    times_diff = np.insert(np.diff(times), 0, 0)

    for i, itime in enumerate(tqdm(times, desc='Determining Orbit and Visit')):
        # if two exposures arent in the same orbit and more than an orbital period apart -> not subsequent orbits but a new visit
        if times_diff[i] * 24 * 60 > 100:
            iorbit_begin = i
            ivisit_begin = i
            iorbit = 0
            ivisit += 1
        # if two exposure more than 10 min apart but less than an orbital period -> subsequent orbits
        elif 10 < times_diff[i] * 24 * 60 <= 100:
            iorbit_begin = i
            iorbit += 1
        # else: two exposures less than 10mins apart -> same orbit and same visit
        iorbits[i] = iorbit
        ivisits[i] = ivisit
        tos[i] = (itime - times[iorbit_begin]) * 24 * 60  # time since first exposure in orbit
        tvs[i] = (itime - times[ivisit_begin]) * 24 * 60  # time since first exposure in visit

    # Only save information of the visits of interest
    mask_visit = [any(tup) for tup in zip(*[ivisits == i for i in meta.which_visits])]
    files = files[mask_visit]
    filter = filter[mask_visit]
    ivisits = ivisits[mask_visit]
    iorbits = iorbits[mask_visit]
    times = times[mask_visit]
    tvs = tvs[mask_visit]
    tos = tos[mask_visit]
    scans = scans[mask_visit]
    exp = exp[mask_visit]

    # changing the numberation of the visits:
    # eg: which_visits = [0,2,5,6]
    # convert ivisits = [0,0,0,0,2,2,2,5,5,5,6,6,6] into [0,0,0,1,1,1,2,2,2,3,3,3]
    ivisits = rankdata(ivisits, method='dense')-1

    print('Writing table into ./filelist.txt')
    table = QTable([files, filter, ivisits, iorbits, times, tvs, tos, scans, exp],
               names=('filenames', 'filter/grism', 'ivisit', 'iorbit', 't_mjd', 't_visit', 't_orbit',
                      'scan', 'exp'))# scan: (0: forward - lower flux, 1: reverse - higher flux, -1: postarg2=0)
    ascii.write(table, meta.workdir + '/filelist.txt', format='ecsv', overwrite=True)

    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    print('Finished s00 \n')

    return meta
