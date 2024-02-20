"""Author: carthik
Revision: 267
Date: 2010-06-08 22:33:22 -0400 (Tue, 08 Jun 2010)
HeadURL: file:///home/esp01/svn/code/python/branches/patricio/photpipe/lib/suntimecorr.py
Id: suntimecorr.py 267 2010-06-09 02:33:22Z carthik
"""
import re
from typing import List
from pathlib import Path

import astropy.time
import numpy as np
from scipy.constants import c

from . import plots
from .options import OPTIONS
from .splinterp import splinterp


def getcoords(file):
    """Use regular expressions to extract X,Y,Z, and time values from the
    horizons file.

    Parameters
    ----------
    file : Strings list
           A list containing the lines of a horizons file.

    Returns
    -------
    A four elements list containing the X, Y, Z, and time arrays of
    values from file.

    Example
    -------
    start_data = '$$SOE'
    end_data   = '$$EOE'

    # Read in whole table as an list of strings, one string per line
    ctable = open('/home/esp01/ancil/horizons/all_spitzer.vec', 'r')
    wholetable = ctable.readlines()
    ctable.close()

    # Find start and end line
    i = 0
    while wholetable[i].find(end_data) == -1:
      if wholetable[i].find(start_data) != -1:
        start = i + 1
      i += 1

    # Chop table
    data = wholetable[start:i-2]

    # Find values:
    x, y, z, t = getcoords(data)

    # print(x, y, z, t)

    Modification History
    --------------------
    2010-11-01  patricio  Written by Patricio Cubillos.
                        pcubillos@fulbrightmail.org
    """
    x, y, z, time = [], [], [], []
    for i in np.arange(len(file)):
        # Use regular expressions to match strings enclosed between X,
        # Y, Z and end of line
        m = re.search(' X =(.*)Y =(.*) Z =(.*)\n', file[i])
        if m is not None:
            x.append(np.double(m.group(1)))
            y.append(np.double(m.group(2)))
            z.append(np.double(m.group(3)))
        # Match first word which is followed by ' = A'
        t = re.search('(.+) = A', file[i])
        if t is not None:
            time.append(np.double(t.group(1)))
    # return numpy arrays
    return np.array(x), np.array(y), np.array(z), np.array(time)


def suntimecorr(meta, obst, coordtable: List[Path], verbose=False):
    """This function calculates the light-travel time correction from
    observer to a standard location.  It uses the 2D coordinates (RA
    and DEC) of the object being observed and the 3D position of the
    observer relative to the standard location.  The latter (and the
    former, for solar-system objects) may be gotten from JPL's
    Horizons system.

    Parameters
    ----------
    meta
                 includes ra, dec and other information
    obst : float or numpy.ndarray
        Time of observation in Julian Date (may be a vector)
    coordtable : str
        Filename of output table from JPL HORIZONS specifying
        the position of the observatory relative to the
        standard position.
    verbose : bool
        If True, print X,Y,Z coordinates.

    Returns
    -------
    This function returns the time correction in seconds to be ADDED
    to the observation time to get the time when the observed photons
    would have reached the plane perpendicular to their travel and
    containing the reference position.

    Notes
    -----
    The position vectors from coordtable are given in the following
    coordinate system:
    Reference epoch: J2000.0
    xy-plane: plane of the Earth's mean equator at the reference epoch
    x-axis  : out along ascending node of instantaneous plane of the Earth's
              orbit and the Earth's mean equator at the reference epoch
    z-axis  : along the Earth mean north pole at the reference epoch

    Ephemerides are often calculated for BJD, barycentric Julian date.
    That is, they are correct for observations taken at the solar
    system barycenter's distance from the target.  The BJD of our
    observation is the time the photons we observe would have crossed
    the sphere centered on the object and containing the barycenter.
    We must thus add the light-travel time from our observatory to
    this sphere.  For non-solar-system observations, we approximate
    the sphere as a plane, and calculate the dot product of the vector
    from the barycenter to the telescope and a unit vector to from the
    barycenter to the target, and divide by the speed of light.

    Properly, the coordinates should point from the standard location
    to the object.  Practically, for objects outside the solar system,
    the adjustment from, e.g., geocentric (RA-DEC) coordinates to
    barycentric coordinates has a negligible effect on the trig
    functions used in the routine.

    The horizons file in coordtable should be in the form of the
    following example, with a subject line of JOB:

    !$$SOF
    !
    ! Example e-mail command file. If mailed to "horizons@ssd.jpl.nasa.gov"
    ! with subject "JOB", results will be mailed back.
    !
    ! This example demonstrates a subset of functions. See main doc for
    ! full explanation. Send blank e-mail with subject "BATCH-LONG" to
    ! horizons@ssd.jpl.nasa.gov for complete example.
    !
     EMAIL_ADDR = 'shl35@cornell.edu'      ! Send output to this address
                                           !  (can be blank for auto-reply)
     COMMAND    = '-79'                  ! Target body, closest apparition

     OBJ_DATA   = 'YES'                    ! No summary of target body data
     MAKE_EPHEM = 'YES'                    ! Make an ephemeris

     START_TIME  = '2005-Aug-24 06:00'     ! Start of table (UTC default)
     STOP_TIME   = '2005-Aug-25 02:00'     ! End of table
     STEP_SIZE   = '1 hour'                 ! Table step-size

     TABLE_TYPE = 'VECTOR'            ! Specify VECTOR ephemeris table type
     CENTER     = '@10'                 ! Set observer (coordinate center)
     REF_PLANE  = 'FRAME'                  ! J2000 equatorial plane

     VECT_TABLE = '3'                      ! Selects output type (3=all).

     OUT_UNITS  = 'KM-S'                   ! Vector units# KM-S, AU-D, KM-D
     CSV_FORMAT = 'NO'                     ! Comma-separated output (YES/NO)
     VEC_LABELS = 'YES'                    ! Label vectors in output (YES/NO)
     VECT_CORR  = 'NONE'                   ! Correct for light-time (LT),
                                           !  or lt + stellar aberration (LT+S),
                                           !  or (NONE) return geometric
                                           !  vectors only.
    !$$EOF

    Example
    -------
    >>> # Spitzer is in nearly the Earth's orbital plane.  Light coming from
    >>> # the north ecliptic pole should hit the observatory and the sun at
    >>> # about the same time.

    >>> import suntimecorr as sc
    >>> ra  = 18.0 * np.pi /  12 # ecliptic north pole coordinates in radians
    >>> dec = 66.5 * np.pi / 180 # "
    >>> obst = np.array([2453607.078])       # Julian date of 2005-08-24 14:00
    >>> print( sc.suntimecorr(ra, dec, obst, '/home/esp01/ancil/horizons/cs41_spitzer.vec') )
    1.00810877 # about 1 sec, close to zero

    >>> # If the object has the RA and DEC of Spitzer, light time should be
    >>> # about 8 minutes to the sun.
    >>> obs  = np.array([111093592.8346969, -97287023.315796047, -42212080.826677799])
    >>> # vector to the object
    >>> obst = np.array([2453602.5])

    >>> print( np.sqrt(np.sum(obs**2.0)) )
    153585191.481 # about 1 AU, good
    >>> raobs  = np.arctan(obs[1]/ obs[0])
    >>> decobs = np.arctan(obs[2]/ np.sqrt(obs[0]**2 + obs[1]**2))
    >>> print(raobs, decobs)
    -0.7192383661, -0.2784282118
    >>> print( sc.suntimecorr(raobs, decobs, obst, '/home/esp01/ancil/horizons/cs41_spitzer.vec') / 60.0)
    8.5228630 # good, about 8 minutes light time to travel 1 AU

    Modification History
    --------------------
    2005-12-01 statia   Written by Statia Luszcz.
    2006-03-09 jh	Corrected 90deg error in algorithm, renamed,
            updated header, made Coordtable a positional
            arg since it's required, switched to radians.
    2007-06-28 jh	Renamed to suntimecorr since we now use
            barycentric Julian date.
    2009-01-28 jh       Change variables to long, use spline instead
            of linfit so we can use one HORIZONS file for
            the whole mission.
    2009-02-22 jh       Reshape spline results to shape of obst.  Make
            it handle unsorted unput data properly.
            Header update.
    2010-07-10 patricio Converted to python. (pcubillos@fulbrightmail.org)
    2010-11-01 patricio Docstring updated.
    December 2021  SZ   Added meta file dependence, check if observations are in vector file and make plots
    """
    ra = meta.ra
    dec = meta.dec

    start_data = '$$SOE'
    end_data = '$$EOE'

    # Read in whole table as an list of strings, one string per line
    with coordtable.open('r', encoding=OPTIONS["encoding"]) as ctable:
        wholetable = ctable.readlines()

    # Find start and end line
    i = 0
    # while end has not been found:
    while wholetable[i].find(end_data) == -1:
        # if start is found get the index of next line:
        if wholetable[i].find(start_data) != -1:
            start = i + 1
        i += 1

    # Chop table
    data = wholetable[start:i - 2]

    # Extract values:
    x, y, z, time = getcoords(data)

    # Interpolate to observing times:
    # We must preserve the shape and order of obst.  Spline takes
    # monotonic input and produces linear output.  x, y, z, time are
    # sorted as HORIZONS produces them.

    # Save shape of obst
    tshape = np.shape(obst)
    #  print(tshape)
    #  print(obst)

    # Reshape to 1D and sort
    obstime = obst.flatten()
    #  print(obstime)

    ########################################## MY PART ###############################

    obs_start = astropy.time.Time(val=obstime[0], format='jd', scale='utc')
    obs_end = astropy.time.Time(val=obstime[-1], format='jd', scale='utc')
    # print('Observation start:', obs_start.iso)
    # print('Observation end:', obs_end.iso)

    horizons_start = astropy.time.Time(val=time[0], format='jd', scale='utc')
    horizons_end = astropy.time.Time(val=time[-1], format='jd', scale='utc')
    # print('Horizon start:', horizons_start.iso)
    # print('Horizon end:', horizons_end.iso)

    if not (min(time) < min(obstime) < max(time)) and (min(time) < max(obstime) < max(time)):
        print('WARNING: HORIZON FILE DOES NOT INCLUDE THE OBSERVATION DATES!!!')
        # exit()
    ########################################## MY PART ###############################
    ti = np.argsort(obstime)  # indexes of sorted array by time
    tsize = np.size(obstime)

    # Allocate output arrays
    obsx = np.zeros(tsize)
    obsy = np.zeros(tsize)
    obsz = np.zeros(tsize)

    # Interpolate sorted arrays
    obsx[ti] = splinterp(obstime[ti], time, x)
    obsy[ti] = splinterp(obstime[ti], time, y)
    obsz[ti] = splinterp(obstime[ti], time, z)

    #  print(obstime[ti], time, x)

    if meta.save_barycorr_plot or meta.show_barycorr_plot:
        plots.barycorr(x, y, z, time, obsx, obsy, obsz, coordtable, meta)

    if verbose:
        print('X, Y, Z = ', obsx, obsy, obsz)

    # Change ra and dec into unit vector n_hat
    object_unit_x = np.cos(dec) * np.cos(ra)
    object_unit_y = np.cos(dec) * np.sin(ra)
    object_unit_z = np.sin(dec)

    # Dot product the vectors with n_hat
    rdotnhat = (obsx * object_unit_x +
                obsy * object_unit_y +
                obsz * object_unit_z)

    # Reshape back to the original shape
    rdotnhat = rdotnhat.reshape(tshape)

    # Time correction is: dt = length/velocity
    # Divide by the speed of light and return
    return rdotnhat / (c / 1000.0)
