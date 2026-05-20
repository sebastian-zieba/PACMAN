import time
from pathlib import Path
from urllib.parse import urlencode
from shutil import copyfileobj
from urllib.request import urlopen
import shutil
import numpy as np
from astropy.io import ascii
from tqdm import tqdm

from .lib import manageevent as me
from .lib import util
from .lib import logedit
from .lib import read_pcf as rd


def run01(pcf_path: Path, meta=None):
    """This function downloads the location of HST during the observations.

    - Retrieves vector data of Hubble from JPL's HORIZONS system on https://ssd.jpl.nasa.gov/horizons_batch.cgi (see Web interface on https://ssd.jpl.nasa.gov/horizons.cgi)
      Based on a perl script found on https://renenyffenegger.ch/notes/Wissenschaft/Astronomie/Ephemeriden/JPL-Horizons
      Also helpful: https://github.com/kevin218/POET/blob/master/code/doc/spitzer_Horizons_README.txt
    - txt file with HST positions in space will be saved in ./run/run_2021-01-01_12-34-56_eventname/ancil/horizons

    .. warning:: This step needs an internet connection!


    Parameters
    ----------
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
        Updated by Sebastian Zieba      May 2025
        New directory structure and logging system
        Modified by Nestor Espinoza     March 2025 
        (with help from ChatGPT)
        Written by Sebastian Zieba      December 2021
    """

    meta, log = util.setup_stage(
        pcf_path=pcf_path,
        stage_num="01",
        previous_stage_num="00",
        copy_filelist=True,
        meta=meta,
    )

    # Read filelist from Stage 00
    filelist_path = meta.workdir / 'filelist.txt'
    filelist = ascii.read(str(filelist_path))

    t_mjd = filelist["t_mjd"]
    ivisit = filelist["ivisit"]

    # Output directory
    horizons_dir = meta.workdir / "ancil" / "horizons"
    horizons_dir.mkdir(parents=True, exist_ok=True)

    base_url = "https://ssd.jpl.nasa.gov/api/horizons.api"

    # Base Horizons parameters. These map closely to the old PACMAN settings.
    base_params = {
        "format": "text",
        "COMMAND": "'-48'",          # Hubble Space Telescope
        "CENTER": "'500@0'",         # Solar System Barycenter
        "MAKE_EPHEM": "'YES'",
        "OBJ_DATA": "'YES'",
        "EPHEM_TYPE": "'VECTORS'",
        "STEP_SIZE": "'5 m'",
        "OUT_UNITS": "'KM-S'",
        "REF_PLANE": "'FRAME'",
        "REF_SYSTEM": "'J2000'",
        "VEC_CORR": "'NONE'",
        "VEC_LABELS": "'YES'",
        "VEC_DELTA_T": "'NO'",
        "CSV_FORMAT": "'NO'",
        "VEC_TABLE": "'3'",
    }

    # save it in ./ancil/bjd_conversion/
    horizons_dir = meta.workdir / 'ancil' / 'horizons'
    horizons_dir.mkdir(parents=True, exist_ok=True)
    for i in tqdm(
        range(int(np.max(ivisit)) + 1),
        desc="Retrieving Horizons file for every visit",
        ascii=True,
    ):
        t_mjd_visit = t_mjd[np.where(ivisit == i)]

        # 1 hour padding before/after the visit, matching the original code
        t_start = np.min(t_mjd_visit) + 2400000.5 - 1.0 / 24.0
        t_end = np.max(t_mjd_visit) + 2400000.5 + 1.0 / 24.0

        params = dict(base_params)
        params["START_TIME"] = f"'JD{t_start:.12f}'"
        params["STOP_TIME"] = f"'JD{t_end:.12f}'"

        # urlencode handles the URL escaping properly.
        url = f"{base_url}?{urlencode(params)}"

        filename = horizons_dir / f"horizons_results_v{i}.txt"

        with urlopen(url) as in_stream, filename.open("wb") as out_file:
            copyfileobj(in_stream, out_file)

    # Save results
    log.writelog('Saving Metadata')
    me.saveevent(meta, meta.workdir / 'WFC3_Meta_Save', save=[])

    log.writelog("Finished s01 \n")
    log.closelog()
    return meta
