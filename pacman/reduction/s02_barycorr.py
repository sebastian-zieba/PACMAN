import os
import numpy as np
from astropy.io import ascii
from astropy.table import Column
from tqdm import tqdm
from ..lib import suntimecorr
from ..lib import util
from ..lib import manageevent as me


def run02(eventlabel, workdir, meta=None):
    """
	Performs the barycentric correction of the observation times

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

	History
	-------
	Written by Sebastian Zieba      December 2021
	"""

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    # read in filelist
    filelist_path = meta.workdir + '/filelist.txt'
    if os.path.exists(filelist_path):
        filelist = ascii.read(filelist_path)

    ivisit = filelist['ivisit']
    t_mjd = filelist['t_mjd']
    t_bjd = np.zeros(len(t_mjd))

    # load in more information into meta
    meta = util.ancil(meta)

    # Converting mjd to bjd
    for i in tqdm(range(max(ivisit) + 1), desc='Converting MJD to BJD'):
        iivisit = ivisit == i
        t_jd = t_mjd[iivisit] + 2400000.5  # converts time to BJD_TDB; see Eastman et al. 2010 equation 4
        t_jd = t_jd + (32.184) / (24.0 * 60.0 * 60.0)
        t_bjd[iivisit] = t_jd + (suntimecorr.suntimecorr(meta, t_jd, meta.coordtable[i], verbose=False)) / (
                60.0 * 60.0 * 24.0)
    # phase = (time - ancil.t0) / ancil.period - math.floor((time - ancil.t0) / ancil.period)
    # if phase > 0.5: phase = phase - 1.0

    print('Writing t_bjd into ./filelist.txt')
    if not any(np.array(filelist.keys()) == 't_bjd'):
        filelist.add_column(Column(data=t_bjd, name='t_bjd'))
        ascii.write(filelist, filelist_path, format='ecsv', overwrite=True)
    else:
        filelist.replace_column(name='t_bjd', col=Column(data=t_bjd, name='t_bjd'))
        ascii.write(filelist, filelist_path, format='rst', overwrite=True)

    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

    print('Finished s02 \n')

    return meta
