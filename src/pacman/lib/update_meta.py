from pathlib import Path

from . import read_pcf as rd
from . import manageevent as me


def update_meta(eventlabel, workdir: Path) -> int:
    """This function reloads the MetaData file from a certain work directory.
    This is needed if a user changed the pcf file in the work directory.

    Notes
    -----
    History:
        Written by Sebastian Zieba      December 2021
    """
    meta = me.loadevent(workdir / f'WFC3_{eventlabel}_Meta_Save')

    pcffile = workdir / 'obs_par.pcf'
    pcf = rd.read_pcf(pcffile)
    rd.store_pcf(meta, pcf)

    meta.workdir = workdir
    me.saveevent(meta, meta.workdir / f'WFC3_{meta.eventlabel}_Meta_Save', save=[])

    print('Successfully reloaded meta file')
    return 0
