# This step creates a table of all available _ima.fits files
# The contents of this filelist.txt are:
# Filename / Filter or Grism / nvisit / norbit / Time since first exposure in orbit / Time since first exposure in visit

import numpy as np
import os, time
import shutil
from astropy.table import QTable
from astropy.io import ascii, fits
from scipy.stats import rankdata
from tqdm import tqdm
from ..lib import read_pcf as rd
from ..lib import util
from ..lib import manageevent as me


class MetaClass:
    """
    A Class which will contain all the metadata of the analysis.
    """
    def __init__(self):
        return


def run00(eventlabel):
    """
    This function does the initial setup of the analysis, including creating a table with information on the observations.

    - 1. Creates a MetaData object
    - 2. Creates a new run directory with the following form, e.g.: ./run/run_2021-01-01_12-34-56_eventname/
    - 3. Copy and pastes the control file (obs_par.pcf) and the fit parameters file (fit_par.txt) into the new directory
    - 4. Reads in all fits files and creates a table which will be saved in filelist.txt.
    - 5. Saves metadata into a file called something like ./run/run_2021-01-01_12-34-56_eventname/WFC3_eventname_Meta_Save.dat

    The information listed in filelist.txt are:

    * **filenames**: The name of the observational file (the file will end with .ima)
    * **instr**: The specific filter or grism used in this observation (taken from the header)
    * **ivisit**: The visit number when the observation was taken (will be calculated in s00)
    * **iorbit**: The orbit number when the observation was taken (will be calculated in s00)
    * **t_mjd**: Mid exposure time (exposure start and end is taken from the header)
    * **t_visit**: Time elapsed since the first exposure in the visit
    * **t_orbit**: Time elapsed since the first exposure in the visit
    * **scan**: Scan direction (0: forward - lower flux, 1: reverse - higher flux, -1: postarg2=0)
    * **exp**: Exposure time

    .. note:: We use the following approach to determine the visit and orbit number:

               - if two exposures arent in the same orbit and more than an orbital period apart -> not subsequent orbits but a new visit
               - if two exposures are more than 10 min apart but less than an orbital period -> subsequent orbits
               - else: two exposures less than 10 mins apart -> same orbit and same visit

    .. note:: We use the following approach to determine the scan direction:

               - if postarg2 < 0   --> scans[i] = 1  --> reverse scan
               - if postarg2 == 0  --> scans[i] = -1 --> no scan direction given
               - else:             --> scans[i] = 0  --> forward scan


    Parameters
    ----------
    eventlabel: str
      the label given to the event in the run script. Will determine the name of the run directory


    Returns
    -------
    meta
      meta object with all the meta data stored in s00


    Notes:
    ----------
    History:
        Written by Sebastian Zieba      December 2021

    """

    print('\nStarting s00')

    # Initialize metadata object
    meta = MetaClass()
    meta.eventlabel = eventlabel

    #this file is saved in /pacman/reduction/s00_table.py
    #topdir is just the path of the directory /pacman/
    meta.pacmandir = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'

    # Create directories for Stage 3 processing
    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
    meta.workdir = 'run_' + datetime + '_' + meta.eventlabel
    if not os.path.exists(meta.workdir):
        os.makedirs(meta.workdir)
    if not os.path.exists(meta.workdir + "/figs"):
        os.makedirs(meta.workdir + "/figs")

    # Load PACMAN control file and store values in Event object
    pcffile = 'obs_par.pcf'
    pcf = rd.read_pcf(pcffile)
    rd.store_pcf(meta, pcf)

    # Copy pcf
    shutil.copy(pcffile, meta.workdir)
    shutil.copy('fit_par.txt', meta.workdir)

    # Create list of file segments
    meta = util.readfiles(meta)
    num_data_files = len(meta.segment_list)
    print(f'Found {num_data_files} data file(s) ending in {meta.suffix}.fits')
    files = meta.segment_list  # gets list of filenames in directory

    # A file called "filelist.txt" should include all _ima.fits file names and if they are a spec or di
    times = np.zeros(len(files))
    exp = np.zeros(len(files))
    instr = np.zeros(len(files), object)
    scans = np.zeros(len(files), dtype=int) #stores scan directions

    # Will create a table with the properties of all _ima.fits files at the first run
    for i, file in enumerate(tqdm(files, desc='Reading in files and their header information')):
        ima = fits.open(file)
        #the header "instr" tells us if the observation used a Filter (-> Direct Image) or a Grism (-> Spectrum)
        instr[i] = ima[0].header['filter']
        exp[i] = ima[0].header['exptime']
        times[i] = (ima[0].header['expstart'] + ima[0].header['expend'])/(2.0)#ima[0].header['expstart']
        #scan direction
        # scan: (0: forward - lower flux, 1: reverse - higher flux, -1: Filter)
        scans[i] = 0  # sets scan direction
        if ima[0].header['postarg2'] < 0: scans[i] = 1
        if instr[i][0] == 'F': scans[i]= -1
        ima.close()

    # files are chronologically sorted
    tsort = np.argsort(times)
    files = np.array([i.split('/')[-1] for i in files])[tsort]
    times = times[tsort]
    exp = exp[tsort]
    instr = instr[tsort]
    instr = np.array([str(iii) for iii in instr])
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

    #TODO: Q: Are there visits with a break of an orbit in between?
    for i, itime in enumerate(tqdm(times, desc='Determining Orbit and Visit')):
        # if two exposures arent in the same orbit and more than an orbital period apart -> not subsequent orbits but a new visit
        if times_diff[i] * 24 * 60 > 100:
            iorbit_begin = i
            ivisit_begin = i
            iorbit = 0
            ivisit += 1
        # if two exposures are more than 10 min apart but less than an orbital period -> subsequent orbits
        elif 10 < times_diff[i] * 24 * 60 <= 100:
            iorbit_begin = i
            iorbit += 1
        # else: two exposures less than 10 mins apart -> same orbit and same visit
        iorbits[i] = iorbit
        ivisits[i] = ivisit
        tos[i] = (itime - times[iorbit_begin]) * 24 * 60  # time since first exposure in orbit
        tvs[i] = (itime - times[ivisit_begin]) * 24 * 60  # time since first exposure in visit


    # TODO: Q: Keep this functionality?
    # Only save information of the visits of interest

    if meta.which_visits != 'everything':
        mask_visit = [any(tup) for tup in zip(*[ivisits == i for i in meta.which_visits])]
        files = files[mask_visit]
        instr = instr[mask_visit]
        ivisits = ivisits[mask_visit]
        iorbits = iorbits[mask_visit]
        times = times[mask_visit]
        tvs = tvs[mask_visit]
        tos = tos[mask_visit]
        scans = scans[mask_visit]
        exp = exp[mask_visit]
        print('The user does not want to analyse every visit (which_visits != everything). '
              'The amount of files analyzed therefore reduced from {0} to {1}.'.format(num_data_files, sum(mask_visit)))


    # changing the numberation of the visits:
    # eg: which_visits = [0,2,5,6]
    # convert ivisits = [0,0,0,0,2,2,2,5,5,5,6,6,6] into [0,0,0,1,1,1,2,2,2,3,3,3]
    ivisits = rankdata(ivisits, method='dense')-1

    print('Writing table into filelist.txt')
    table = QTable([files, instr, ivisits, iorbits, times, tvs, tos, scans, exp],
               names=('filenames', 'instr', 'ivisit', 'iorbit', 't_mjd', 't_visit', 't_orbit',
                      'scan', 'exp'))# scan: (0: forward - lower flux, 1: reverse - higher flux, -1: postarg2=0)
    ascii.write(table, meta.workdir + '/filelist.txt', format='rst', overwrite=True)

    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    print('Finished s00 \n')

    return meta
