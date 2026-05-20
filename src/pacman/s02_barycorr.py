import time
from pathlib import Path
import shutil
import numpy as np

from astropy.io import ascii
from astropy.table import Column
from tqdm import tqdm

from .lib import suntimecorr
from .lib import util
from .lib import manageevent as me
from .lib import logedit


def run02(pcf_path: Path, meta=None):
    """Performs the barycentric correction of the observation times

    - performs the barycentric correction based on the t_mjd in filelist.txt.
    - Adds another column to filelist.txt called t_bjd
    - Plots will be saved in ./run/run_2021-01-01_12-34-56_eventname/ancil/horizons

    Parameters
    ----------
    eventlabel : str
        the label given to the event in the run script. Will determine the name of the run directory
    workdir : str
        the name of the work directory.
    meta
        the name of the metadata file

    Returns
    -------
    meta
        meta object with all the meta data stored in s01

    Notes
    -----
    History:
        Written by Sebastian Zieba      December 2021
    """

    pcf_path = Path(pcf_path)
    rundir = pcf_path.parent

    # Find latest Stage 01 workdir
    s01_workdir = util.find_latest_stage_run(rundir, 'stage01', 's01_run_*')

    if meta is None:
        meta = me.loadevent(s01_workdir / 'WFC3_Meta_Save')

    # Create new Stage 02 workdir
    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
    meta.inputdir = s01_workdir
    meta.stage02dir = rundir / 'stage02'
    meta.workdir = meta.stage02dir / f's02_run_{datetime}'
    meta.workdir.mkdir(parents=True, exist_ok=True)

    # Copy current config files from pacman_run_files
    shutil.copy(pcf_path / 'obs_par.pcf', meta.workdir)
    shutil.copy(pcf_path / 'fit_par.txt', meta.workdir)

    # Copy filelist from previous stage
    shutil.copy(meta.inputdir / 'filelist.txt', meta.workdir)

    previous_log = meta.inputdir / "s01.log"
    meta.logname = meta.workdir / "s02.log"
    log = logedit.Logedit(meta.logname, read=previous_log)

    log.writelog("Starting s02")
    log.writelog(f"Using Stage 01 input directory: {meta.inputdir}")
    log.writelog(f"Location of the new Stage 02 run directory: {meta.workdir}")

    # Read filelist from Stage 00/01 input lineage
    filelist_path = meta.workdir / 'filelist.txt'
    filelist = ascii.read(str(filelist_path))

    ivisit = filelist['ivisit']
    t_mjd = filelist['t_mjd']
    t_bjd = np.zeros(len(t_mjd))

    # load in more information into meta
    meta = util.ancil(meta)

    # Converting mjd to bjd
    for i in tqdm(range(max(ivisit) + 1), desc='Converting MJD to BJD', ascii=True):
        iivisit = ivisit == i
        t_jd = t_mjd[iivisit] + 2400000.5  # converts time to BJD_TDB; see Eastman et al. 2010 equation 4
        t_jd = t_jd + (32.184) / (24.0 * 60.0 * 60.0)
        t_bjd[iivisit] = t_jd + (suntimecorr.suntimecorr(meta, t_jd, meta.coordtable[i], verbose=False)) / (
                60.0 * 60.0 * 24.0)

    # Identify orbits and visits
    iorbit = 0
    ivisit = 0
    tos = np.zeros(len(t_bjd))  # Time since begin of orbit
    tvs = np.zeros(len(t_bjd))  # Time since begin of visit
    iorbit_begin = 0  # Index of first exposure in orbit
    ivisit_begin = 0  # Index of first exposure in visit
    times_diff = np.insert(np.diff(t_bjd), 0, 0)

    for i, itime in enumerate(tqdm(t_bjd, desc='Correcting orbit and visit times to BJD', ascii=True)):
        # If two exposures arent in the same orbit and more than an orbital period apart -> not subsequent orbits but a new visit
        if times_diff[i] * 24 * 60 > 100:
            iorbit_begin = i
            ivisit_begin = i
            iorbit = 0
            ivisit += 1
        # If two exposures are more than 10 min apart but less than an orbital period -> subsequent orbits
        elif 30 < times_diff[i] * 24 * 60 <= 100:
            iorbit_begin = i
            iorbit += 1
        # Else: two exposures less than 10 mins apart -> same orbit and same visit
        tos[i] = (itime - t_bjd[iorbit_begin]) * 24 * 60  # time since first exposure in orbit
        tvs[i] = (itime - t_bjd[ivisit_begin]) * 24 * 60  # time since first exposure in visit

    log.writelog('Writing t_bjd into filelist.txt')

    if 't_bjd' not in filelist.colnames:
        filelist.add_column(Column(data=t_bjd, name='t_bjd'))
    else:
        filelist.replace_column(name='t_bjd', col=Column(data=t_bjd, name='t_bjd'))

    # Overwrite old visit and orbit times with BJD corrected ones
    filelist['t_visit'] = tvs
    filelist['t_orbit'] = tos

    # Save updated filelist into the Stage 02 workdir
    ascii.write(filelist, meta.workdir / 'filelist.txt', format='rst', overwrite=True)

    log.writelog('Saving Metadata')
    me.saveevent(meta, meta.workdir / 'WFC3_Meta_Save', save=[])

    log.writelog('Finished s02 \n')
    log.closelog()

    return meta