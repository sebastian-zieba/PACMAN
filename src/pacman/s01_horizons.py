from pathlib import Path
from shutil import copyfileobj
from urllib.request import urlopen

import numpy as np
from astropy.io import ascii
from tqdm import tqdm

from .lib import manageevent as me
from .lib.options import OPTIONS


def run01(eventlabel, workdir: Path, meta=None):
    """This function downloads the location of HST during the observations.

    - Retrieves vector data of Hubble from JPL's HORIZONS system on https://ssd.jpl.nasa.gov/horizons_batch.cgi (see Web interface on https://ssd.jpl.nasa.gov/horizons.cgi)
      Based on a perl script found on https://renenyffenegger.ch/notes/Wissenschaft/Astronomie/Ephemeriden/JPL-Horizons
      Also helpful: https://github.com/kevin218/POET/blob/master/code/doc/spitzer_Horizons_README.txt
    - txt file with HST positions in space will be saved in ./run/run_2021-01-01_12-34-56_eventname/ancil/horizons

    .. warning:: This step needs an internet connection!


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
       meta object with all the meta data stored in s00

    Notes:
    ----------
    History:
        Written by Sebastian Zieba      December 2021
    """
    print('Starting s01')
    if meta is None:
        meta = me.loadevent(workdir / f'WFC3_{eventlabel}_Meta_Save')

    # read in filelist
    filelist_path = meta.workdir / 'filelist.txt'
    if filelist_path.exists():
        filelist = ascii.read(filelist_path)

    t_mjd = filelist['t_mjd']
    ivisit = filelist['ivisit']

    # Setting for the JPL Horizons interface
    settings = [
        "COMMAND= -48",  # Hubble
        "CENTER= 500@0",  # Solar System Barycenter (SSB) [500@0]
        "MAKE_EPHEM= YES",
        "TABLE_TYPE= VECTORS",
        # "START_TIME= $ARGV[0]",
        # "STOP_TIME= $ARGV[1]",
        "STEP_SIZE= 5m",  # 5 Minute interval
        "OUT_UNITS= KM-S",
        "REF_PLANE= FRAME",
        "REF_SYSTEM= J2000",
        "VECT_CORR= NONE",
        "VEC_LABELS= YES",
        "VEC_DELTA_T= NO",
        "CSV_FORMAT= NO",
        "OBJ_DATA= YES",
        "VEC_TABLE= 3"]

    # Replacing symbols for URL encoding
    for i, _ in enumerate(settings):
        settings[i] = settings[i].replace(" =", "=").replace("= ", "=")
        settings[i] = settings[i].replace(" ", "%20")
        settings[i] = settings[i].replace("&", "%26")
        settings[i] = settings[i].replace(";", "%3B")
        settings[i] = settings[i].replace("?", "%3F")

    settings = '&'.join(settings)
    settings = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&' + settings

    # save it in ./ancil/bjd_conversion/
    horizons_dir = meta.workdir / 'ancil' / 'horizons'
    if not horizons_dir.exists():
        horizons_dir.mkdir(parents=True)

    # retrieve positions for every individual visit
    for i in tqdm(range(max(ivisit) + 1), desc='Retrieving Horizons file for every visit', ascii=True):
        t_mjd_visit = t_mjd[np.where(ivisit == i)]
        t_start = min(
            t_mjd_visit) + 2400000.5 - 1 / 24  # Start of Horizons file one hour before first exposure in visit
        t_end = max(t_mjd_visit) + 2400000.5 + 1 / 24  # End of Horizons file one hour after last exposure in visit

        # Complete the settings by also adding information on the start and end of the times of interest.
        set_start = "START_TIME=JD{0}".format(t_start)
        set_end = "STOP_TIME=JD{0}".format(t_end)

        # Full link
        settings_new = settings + '&' + set_start + '&' + set_end

        # Location where to save the data
        filename = horizons_dir / f'horizons_results_v{i}.txt'

        # Download data
        with urlopen(settings_new) as in_stream, filename.open('wb') as out_file:
            copyfileobj(in_stream, out_file)

    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir / f'WFC3_{meta.eventlabel}_Meta_Save', save=[])

    print('Finished s01 \n')
    return meta
