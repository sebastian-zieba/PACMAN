import os
import numpy as np
from astropy.io import ascii
from astropy.table import Column
from tqdm import tqdm
from ..lib import suntimecorr
from ..lib import util
from ..lib import manageevent as me


def run02(eventlabel, workdir, meta=None):

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
	for i in tqdm(range(max(ivisit)+1), desc='Converting MJD to BJD'):
		iivisit = ivisit == i
		t_jd = t_mjd[iivisit] + 2400000.5  # converts time to BJD_TDB; see Eastman et al. 2010 equation 4
		t_jd = t_jd + (32.184) / (24.0 * 60.0 * 60.0)
		t_bjd[iivisit] = t_jd + (suntimecorr.suntimecorr(meta, t_jd, meta.coordtable[i], verbose=False)) / (60.0 * 60.0 * 24.0)
		#phase = (time - ancil.t0) / ancil.period - math.floor((time - ancil.t0) / ancil.period)
		#if phase > 0.5: phase = phase - 1.0

	print('Writing t_bjd into ./filelist.txt')
	if not any(np.array(filelist.keys())=='t_bjd'):
		filelist.add_column(Column(data = t_bjd, name='t_bjd'))
		ascii.write(filelist, filelist_path, format='ecsv', overwrite=True)
	else:
		filelist.replace_column(name='t_bjd', col=Column(data = t_bjd, name='t_bjd'))
		ascii.write(filelist, filelist_path, format='ecsv', overwrite=True)

	# Save results
	print('Saving Metadata')
	me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

	print('Finished s02 \n')

	return meta
