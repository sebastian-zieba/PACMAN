# This code computes the mean position of the direct image for each visit
import sys
sys.path.insert(0, './util')

from pylab import *
import gaussfitter, ancil, os, glob
from astropy.io import fits
import yaml
from astropy.io import ascii

yaml_path = "config/obs_par.yaml"
with open(yaml_path, 'r') as file:
	params = yaml.safe_load(file)

ancil = ancil.AncillaryData(params)

if ancil.output ==True:
	f = open('config/xrefyref.txt', 'w')						#opens file to store positions of reference pixels

filetable = ascii.read('config/filelist.txt')

mask_di = filetable['filter/grism'] == ancil.FILTER

files_di = [ancil.path + '/' + i for i in filetable['filenames'][mask_di].data]
print(files_di)

for i in enumerate(files_di):
	ima = fits.open(i[1])
	print(ima[0].header['OBSTYPE'])

	if (ima[0].header['filter'] == ancil.FILTER):
		print("filename", [i[1]])
		LTV1 = ima[1].header['LTV1']					#X offset to get into physical pixels
		LTV2 = ima[1].header['LTV2']					#Y offset to get to physical pixels
		nrow = len(ima[1].data[:,0])
		ncol = len(ima[1].data[0,:])
		t = ima[0].header['expstart']

		dat = ima[1].data[ancil.rmin:ancil.rmax, ancil.cmin:ancil.cmax]				#cuts out stamp around the target star
		err = ima[2].data[ancil.rmin:ancil.rmax, ancil.cmin:ancil.cmax]

		results = gaussfitter.gaussfit(dat, err)
			
		if ancil.diagnostics==True:
			plt.title("Direct image")
			plt.imshow(dat*ima[0].header['exptime'], origin ='lower',vmin=0, vmax=5000)
			plt.plot(results[2], results[3],marker='x', color='orange', markeredgewidth=3., ms=10, label='centroid', linestyle="none")
			plt.legend(numpoints=1)
			plt.colorbar()
			if ancil.output == True:
				if not os.path.isdir('config/images/'):
					os.makedirs('config/images/')
				plt.savefig('config/images/{0}.png'.format(i[0]))
			plt.show()
		
		if ancil.output==True:
			print(t, results[3]+ancil.rmin-LTV1, results[2]+ancil.cmin-LTV2, file=f)

		print(t, results[3]+ancil.rmin-LTV1, results[2]+ancil.cmin-LTV2, [i[1]]) 		#fit results
	ima.close()

if ancil.output==True:
	f.close()
