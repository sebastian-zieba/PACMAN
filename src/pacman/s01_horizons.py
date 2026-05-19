from pathlib import Path
from urllib.parse import urlencode
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
        Modified by Nestor Espinoza     March 2025 
        (with help from ChatGPT)
        Written by Sebastian Zieba      December 2021
    """

    print("Starting s01")

    if meta is None:
        meta = me.loadevent(workdir / f"WFC3_{eventlabel}_Meta_Save")

    # Read filelist
    filelist_path = meta.workdir / "filelist.txt"
    if not filelist_path.exists():
        raise FileNotFoundError(f"Could not find {filelist_path}")

    filelist = ascii.read(filelist_path)
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
        params["START_TIME"] = f"'JD{t_start:.8f}'"
        params["STOP_TIME"] = f"'JD{t_end:.8f}'"

        # urlencode handles the URL escaping properly.
        url = f"{base_url}?{urlencode(params)}"

        filename = horizons_dir / f"horizons_results_v{i}.txt"

        with urlopen(url) as in_stream, filename.open("wb") as out_file:
            copyfileobj(in_stream, out_file)

    print("Saving Metadata")
    me.saveevent(meta, meta.workdir / f"WFC3_{meta.eventlabel}_Meta_Save", save=[])

    print("Finished s01\n")
    return meta
