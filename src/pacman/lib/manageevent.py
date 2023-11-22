"""Name
  ----
  Manage Event

  File
  ----
  manageevnet.py

  Description
  -----------
  Routines for handling events.

  Package Contents
  ----------------
  saveevent(event, filename, save=['event'], delete=[])
      Saves an event in .dat (using cpickle) and .h5 (using h5py) files.

  loadevent(filename, load):
      Loads an event stored in .dat and .h5 files.

  updateevent(event, filename, add):
      Adds parameters given by add from filename to event.


  Examples
  --------
  >>> from manageevent import *

  # Save  hd209bs51_ini.dat and hd209bs51_ini.h5 files.
  >>> saveevent(event, 'd209bs51_ini',
                save=['data', 'head','uncd', 'bdmskd'])

  # Load the event and its data frames
  >>> event = loadevent('hd209bs51_ini', ['data'])

  # Load uncd and bdmsk into event:
  >>> updateevent(event, 'hd209bs51_ini', ['uncd', 'bdmskd'])

  Revisions
  ---------
  2010-07-10  patricio  joined loadevent and       pcubillos@fulbrightmail.org
                        saveevent into this package.
                        updateevent added.
  2010-11-12  patricio  reimplemented using exec()
"""
from pathlib import Path
from typing import List, Optional

import h5py as h5

from .options import OPTIONS

try:
    import cPickle as pickle
except Exception:
    import _pickle as pickle


def saveevent(event, filename: Path, save: Optional[List[str]] = [],
              delete: Optional[List[str]] = [], protocol: Optional[int] = 3):
    """Saves an event in .dat (using cpickle) and .h5 (using h5py) files.

    Parameters
    ----------
    event :
        An Event instance.
    filename : pathlib.Path
        Path to the event file.
    save : list of str, optional
        The elements of this tuple contain the parameters to save.
        We usually use the values: 'data', 'uncd', 'head', 'bdmskd',
        'brmksd' or 'mask'.
    delete : list of str, optional
        Parameters to be deleted.

    Notes
    -----
    The input filename should not have the .dat nor the .h5 extentions.
    Side effect: This routine deletes all parameters except 'event' after saving it.

    Examples
    --------
    See package example.

    Revisions
    ---------
    2010-07-10  patricio  Added documentation.     pcubillos@fulbrightmail.org
    """
    filename = filename if isinstance(filename, Path) else filename
    if save != []:
        with h5.File(f'{filename}.h5', 'w') as handle:
            for param in save:
                handle[param] = getattr(event, param)
                delattr(event, param)

                # NOTE: Calibration data
                if event.havecalaor:
                    handle[f"pre{param}"] = getattr(event, f"pre{param}")
                    handle[f"post{param}"] = getattr(event, f"post{param}")
                    delattr(event, f"pre{param}")
                    delattr(event, f"post{param}")

    # NOTE: Delete if requested
    for param in delete:
        delattr(event, param)
        if event.havecalaor:
            delattr(event, f"pre{param}")
            delattr(event, f"post{param}")

    # NOTE: Pickle-Save the event
    with Path(f'{filename}.dat').open('wb') as handle:
        pickle.dump(event, handle, protocol)


def loadevent(filename: Path, load: Optional[List[str]] = [],
              loadfilename: Optional[bool] = None):
    """Loads an event stored in .dat and .h5 files.

    Parameters
    ----------
    filename : pathlib.Path
        Path to the event file.
    load : list of str
        The elements of this tuple contain the parameters to read.
        We usually use the values: 'data', 'uncd', 'head', 'bdmskd',
        'brmskd' or 'mask'.

    Notes
    -----
    The input filename should not have the .dat nor the .h5 extentions.

    Returns
    -------
    This function return an Event instance.

    Examples
    --------
    See package example.

    Revisions
    ---------
    2010-07-10  patricio  Added documentation.     pcubillos@fulbrightmail.org
    """
    filename = filename if isinstance(filename, Path) else filename
    with Path(f'{filename}.dat').open('rb') as handle:
        event = pickle.load(handle, encoding=OPTIONS["encoding"])

    if loadfilename is None:
        loadfilename = filename

    if load:
        with h5.File(f'{loadfilename}.h5', 'r') as handle:
            for param in load:
                setattr(event, param, handle[param][:])

                # NOTE: Calibration data
                if event.havecalaor:
                    setattr(even, f"pre{param}", handle[f"pre{param}"][:])
                    setattr(even, f"post{param}", handle[f"post{param}"][:])
    return event


def updateevent(event, filename: Path, add: List[str]):
    """Adds parameters given by add from filename to event.

    Parameters
    ----------
    event : An Event instance.
    filename : pathlib.Path
        Path to the event file.
    add : list of str
        The elements of this tuple contain the parameters to
        add.  We usually use the values: 'data', 'uncd', 'head',
        'bdmskd', 'brmaskd' or 'mask'.

    Notes
    -----
    The input filename should not have the .dat nor the .h5 extentions.

    Returns
    -------
    This function return an Event instance.

    Examples
    --------
    See package example.

    Revisions
    ---------
    2010-07-10  patricio  Created.     pcubillos@fulbrightmail.org
    """
    event2 = loadevent(filename, load=add)
    for param in add:
        setattr(event, param, getattr(event2, param))

        # NOTE: Calibration data
        if event.havecalaor:
            setattr(event, f"pre{param}", getattr(event2, f"pre{param}"))
            setattr(event, f"post{param}", getattr(event2, f"post{param}"))
    return event
