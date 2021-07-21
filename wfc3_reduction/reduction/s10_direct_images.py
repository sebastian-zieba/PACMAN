# This code computes the mean position of the direct image for each visit
import sys

from ..lib import gaussfitter

import matplotlib.pyplot as plt
import os
from ..lib import manageevent as me
from astropy.io import ascii, fits

from tqdm import tqdm


def run10(eventlabel, workdir, meta=None):

	if meta == None:
		meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

	if meta.direct_image_output == True:
		f = open(meta.workdir + '/xrefyref.txt', 'w')						#opens file to store positions of reference pixels

	filelist_path = ascii.read(meta.workdir + '/filelist.txt')

	mask_di = filelist_path['filter/grism'] == meta.filter
	t_bjd = filelist_path['t_bjd'][mask_di]
	files_di = [meta.path + '/' + i for i in filelist_path['filenames'][mask_di].data]


	#iterate over the direct images
	for i, file in enumerate(files_di):

		ima = fits.open(file)
		#print(ima[0].header['OBSTYPE'])

		if (ima[0].header['filter'] == meta.filter):
			#print("filename", [file])
			LTV1 = ima[1].header['LTV1']					#X offset to get into physical pixels
			LTV2 = ima[1].header['LTV2']					#Y offset to get to physical pixels
			nrow = len(ima[1].data[:,0])
			ncol = len(ima[1].data[0,:])
			t = t_bjd[i]
			#print(ima[1].data.shape)

			dat = ima[1].data[meta.di_rmin:meta.di_rmax, meta.di_cmin:meta.di_cmax]				#cuts out stamp around the target star
			err = ima[2].data[meta.di_rmin:meta.di_rmax, meta.di_cmin:meta.di_cmax]

			results = gaussfitter.gaussfit(dat, err)

			if meta.direct_image_diagnostics==True:
				plt.clf()
				plt.title("Direct image #{0}, nvisit {1}, norbit {2}".format(i, filelist_path['nvisit'][mask_di][i], filelist_path['norbit'][mask_di][i]))
				plt.imshow(dat*ima[0].header['exptime'], origin ='lower',vmin=0, vmax=5000)
				plt.plot(results[2], results[3],marker='x', color='orange', markeredgewidth=3., ms=10, label='centroid', linestyle="none")
				plt.legend(numpoints=1)
				plt.colorbar()
				if meta.direct_image_output == True:
					if not os.path.isdir(meta.workdir + '/figs/images/'):
						os.makedirs(meta.workdir + '/figs/images/')
					plt.savefig(meta.workdir + '/figs/images/{0}.png'.format(i))
				plt.show()
				plt.close()

			if meta.direct_image_output==True:
				print(t, results[3]+meta.di_rmin-LTV1, results[2]+meta.di_cmin-LTV2, file=f)

			#print(t, results[3]+ancil.rmin-LTV1, results[2]+ancil.cmin-LTV2, [file]) 		#fit results
		ima.close()

	if meta.direct_image_output==True:
		f.close()

	print('Finished 10.py')

	return meta
