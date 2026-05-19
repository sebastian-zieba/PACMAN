import time
from importlib import resources
import shutil
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.table import QTable
from astropy.io import ascii, fits
from scipy.stats import rankdata
from tqdm import tqdm

from .lib import read_pcf as rd
from .lib import util
from .lib import manageevent as me
from .lib import plots


class MetaClass:
    """A Class which will contain all the metadata of the analysis."""
    def __init__(self):
        return


def run00(pcf_path: Optional[Path] = Path.cwd()):
    """This function does the initial setup of the analysis, including
    creating a table with information on the observations. This table will
    be saved into 'filelist.txt'.

    Steps:

    - 1. Creates a metadata object.
    - 2. Reads the live ``obs_par.pcf`` file from ``pacman_run_files``.
    - 3. Creates a new Stage 00 run directory with the form
         ``rundir/stage00/s00_run_YYYY-MM-DD_HH-MM-SS``.
    - 4. Copies ``obs_par.pcf`` and ``fit_par.txt`` into the Stage 00 run directory
         as provenance snapshots.
    - 5. Reads all FITS files and creates ``filelist.txt``.
    - 6. Saves metadata into the Stage 00 run directory.

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


    Returns
    -------
    meta
        Meta object with all the meta data stored in s00


    Notes:
    ----------
    History:
        Written by Sebastian Zieba      December 2021
    """
    print('\nStarting s00')

    # Initialize metadata object
    meta = MetaClass()

    # The live/current PACMAN control files are always read from pacman_run_files.
    # Example:
    #   rundir/pacman_run_files/obs_par.pcf
    #   rundir/pacman_run_files/fit_par.txt
    pcf_path = Path(pcf_path)
    pcffile = pcf_path / 'obs_par.pcf'
    fit_parfile = pcf_path / 'fit_par.txt'

    # Load PACMAN control file and store values in metadata object
    pcf = rd.read_pcf(pcffile)
    rd.store_pcf(meta, pcf)

    # PACMAN package directory
    meta.pacmandir = resources.files("pacman")
    print('Location of PACMAN:', meta.pacmandir)

    # Create the Stage 00 output directory:
    #   rundir/stage00/s00_run_YYYY-MM-DD_HH-MM-SS
    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')

    meta.stage00dir = meta.rundir / 'stage00'
    meta.workdir = meta.stage00dir / f's00_run_{datetime}'

    meta.workdir.mkdir(parents=True, exist_ok=True)
    print('Location of the new Stage 00 run directory:', meta.workdir)

    # Create figure directory inside this Stage 00 run
    figure_dir = meta.workdir / 'figs'
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Copy obs_par.pcf and fit_par.txt into the Stage 00 run directory as provenance snapshots.
    # These copied files document the settings used for this run.
    # Future runs still read the live files from pacman_run_files.
    shutil.copy(pcffile, meta.workdir)
    shutil.copy(fit_parfile, meta.workdir)
    print('pcf and fit_par files copied to the Stage 00 run directory:', meta.workdir)

    # Create list of file segments
    meta = util.readfiles(meta)
    num_data_files = len(meta.segment_list)
    print(f'Found {num_data_files} data file(s) ending in {meta.suffix}.fits')
    files = meta.segment_list  # gets list of filenames in directory

    # A file called "filelist.txt" should include all _ima.fits file names and if they are a spec or di
    times = np.zeros(len(files))
    exp = np.zeros(len(files))
    instr = np.zeros(len(files), object)
    scans = np.zeros(len(files), dtype=int)  # Stores scan directions

    postarg2 = np.zeros(len(files))

    # Will create a table with the properties of all _ima.fits files at the first run
    for i, file in enumerate(tqdm(files, desc='Reading in files and their headers', ascii=True)):
        ima = fits.open(file)
        # The header "instr" tells us if the observation used a Filter (-> Direct Image) or a Grism (-> Spectrum)
        # print(repr(ima[0].header))
        instr[i] = ima[0].header['filter']
        exp[i] = ima[0].header['exptime']
        times[i] = (ima[0].header['expstart'] + ima[0].header['expend'])/(2.0) # mid exposure time
        # Scan direction
        # Scan: (0: forward - lower flux, 1: reverse - higher flux, -1: Filter)
        scans[i] = 0  # sets scan direction
        if ima[0].header['postarg2'] < 0: scans[i] = 1
        if instr[i][0] == 'F': scans[i]= -1
        postarg2[i] = ima[0].header['postarg2']
        ima.close()

    # Files are chronologically sorted
    tsort = np.argsort(times)
    files = np.array([Path(file).name for file in files])[tsort]
    # files = np.array(files)
    # files = files[tsort]
    times = times[tsort]
    exp = exp[tsort]
    instr = instr[tsort]
    instr = np.array([str(iii) for iii in instr])
    scans = scans[tsort]
    postarg2 = postarg2[tsort]
    #
    # for i, file in enumerate(tqdm(files, desc='Reading in files and their headers', ascii=True)):
    #     f = open(meta.workdir + '/{0}header_{1}.txt'.format(i, str(file).split('/')[-1]), 'w')
    #     #print(file)
    #     ima = fits.open(file)
    #     print(repr(ima[0].header), file=f)
    #     ima.close()
    #     f.close

    # Identify orbits and visits
    iexp_orb = np.zeros(len(times), dtype=int)
    iorbits = np.zeros(len(times), dtype=int)
    ivisits = np.zeros(len(times), dtype=int)

    iorbit = 0
    ivisit = 0
    tos = np.zeros(len(times)) #time since begin of orbit
    tvs = np.zeros(len(times)) #time since begin of visit
    iexp_orb_counter = 0 #index of first exposure in orbit
    iorbit_begin = 0 #index of first exposure in orbit
    ivisit_begin = 0 #index of first exposure in visit
    times_diff = np.insert(np.diff(times), 0, 0)

    for i, itime in enumerate(tqdm(times, desc='Determining orbit(s) and visit(s)', ascii=True)):
        # if two exposures arent in the same orbit and more than an orbital period apart -> not subsequent orbits but a new visit
        if times_diff[i] * 24 * 60 > 100:
            iorbit_begin = i
            ivisit_begin = i
            iorbit = 0
            ivisit += 1
            iexp_orb_counter = 0
        # if two exposures are more than 10 min apart but less than an orbital period -> subsequent orbits
        elif 30 < times_diff[i] * 24 * 60 <= 100:
            iorbit_begin = i
            iorbit += 1
            iexp_orb_counter = 0
        # else: two exposures less than 10 mins apart -> same orbit and same visit
        iorbits[i] = iorbit
        ivisits[i] = ivisit
        tos[i] = (itime - times[iorbit_begin]) * 24 * 60  # time since first exposure in orbit
        tvs[i] = (itime - times[ivisit_begin]) * 24 * 60  # time since first exposure in visit
        if instr[i][0] == 'F':
            iexp_orb[i] = -1
            iexp_orb_counter = 0
        elif instr[i][0] == 'G':
            iexp_orb[i] = iexp_orb_counter
            iexp_orb_counter += 1

    # Saves a plot showing when was observered
    if meta.save_obs_times_plot or meta.show_obs_times_plot:
        plots.obs_times(meta, times, ivisits, iorbits)

    # Only save information of the visits of interest
    if meta.which_visits != 'everything':
        mask_visit = [any(tup) for tup in zip(*[ivisits == i for i in meta.which_visits])]
        files = files[mask_visit]
        instr = instr[mask_visit]
        iexp_orb = iexp_orb[mask_visit]
        ivisits = ivisits[mask_visit]
        iorbits = iorbits[mask_visit]
        times = times[mask_visit]
        tvs = tvs[mask_visit]
        tos = tos[mask_visit]
        scans = scans[mask_visit]
        postarg2 = postarg2[mask_visit]
        exp = exp[mask_visit]
        print('The user does not want to analyse every visit (which_visits != everything). '
              'The amount of files analyzed therefore reduced from {0} to {1}.'.format(num_data_files, sum(mask_visit)))

    # changing the numeration of the visits:
    # eg: when which_visits = [0,2,5,6]
    # then the ivisits list will look something like this:
    # ivisits = [0,0,0,0,2,2,2,5,5,5,6,6,6]
    # the following will then reindex ivisit into [0,0,0,1,1,1,2,2,2,3,3,3]
    ivisits = rankdata(ivisits, method='dense')-1

    # Saves a plot showing when was observed with the new indexing of ivisits if which_visits was set by the user
    if meta.which_visits  != 'everything' and (meta.save_obs_times_plot or meta.show_obs_times_plot):
        plots.obs_times(meta, times, ivisits, iorbits, updated=True)

    print('Writing table into filelist.txt')
#    table = QTable([files, instr, ivisits, iorbits, iexp_orb, times, tvs, tos, scans, exp, postarg2],
#               names=('filenames', 'instr', 'ivisit', 'iorbit', 'iexp_orb', 't_mjd', 't_visit', 't_orbit',
#                      'scan', 'exp', 'postarg2'))# scan: (0: forward - lower flux, 1: reverse - higher flux, -1: Direct Image)
    table = QTable([files, instr, ivisits, iorbits, iexp_orb, times, tvs, tos, scans, exp],
               names=('filenames', 'instr', 'ivisit', 'iorbit', 'iexp_orb', 't_mjd', 't_visit', 't_orbit',
                      'scan', 'exp'))# scan: (0: forward - lower flux, 1: reverse - higher flux, -1: Direct Image)
    ascii.write(table, meta.workdir / 'filelist.txt', format='rst', overwrite=True)

    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir / 'WFC3_Meta_Save', save=[])

    print('Finished s00 \n')
    return meta
