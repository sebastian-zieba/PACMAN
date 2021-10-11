# This code computes the mean position of the direct image for each visit

from ..lib import gaussfitter
from ..lib import manageevent as me
from astropy.io import ascii, fits
from ..lib import plots
from tqdm import tqdm
import os
from ..lib import util


def run10(eventlabel, workdir, meta=None):

	if meta == None:
		meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

	f = open(meta.workdir + '/xrefyref.txt', 'w')						#opens file to store positions of reference pixels

	# read in filelist
	filelist_path = meta.workdir + '/filelist.txt'
	if os.path.exists(filelist_path):
		filelist = ascii.read(filelist_path)

	# load in more information into meta
	meta = util.ancil(meta, s10=True)


	t_bjd = meta.t_bjd_di


	#iterate over the direct images
	for i, file in enumerate(tqdm(meta.files_di, desc='Determining Source Positions for Direct Images')):

		ima = fits.open(file)

		LTV1 = ima[1].header['LTV1']					#X offset to get into physical pixels
		LTV2 = ima[1].header['LTV2']					#Y offset to get to physical pixels

		t = t_bjd[i]

		dat = ima[1].data[meta.di_rmin:meta.di_rmax, meta.di_cmin:meta.di_cmax]				#cuts out stamp around the target star
		err = ima[2].data[meta.di_rmin:meta.di_rmax, meta.di_cmin:meta.di_cmax]

		plots.image_quick(ima, i, meta)

		results = gaussfitter.gaussfit(dat, err)

		if meta.save_image_plot or meta.show_image_plot:
			plots.image(dat, ima, results, i, meta)

		print(t, results[3]+meta.di_rmin-LTV1, results[2]+meta.di_cmin-LTV2, file=f)

		ima.close()

	f.close()

	# Save results
	print('Saving Metadata')
	me.saveevent(meta, meta.workdir + '/WFC3_' + meta.eventlabel + "_Meta_Save", save=[])

	print('Finished s10 \n')

	return meta
