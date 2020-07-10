# This code computes the mean position of the direct image for each visit
import sys
sys.path.insert(0, './util')
from pylab import *
import gaussfitter, pyfits, os, glob, math
from astropy.io import ascii

def convert_to_bool(str):
	if str=="True": return True
	elif str=="False": return False
	else: return "String not equal to True or False"

def make_dict(table):
	return {x['parameter']: x['value'] for x in table}
	
obs_par = make_dict(ascii.read("config/obs_par.txt", Reader=ascii.CommentedHeader))
path = obs_par['path'] 

output = convert_to_bool(obs_par['direct_image_output']) 			#flag specifying whether coordinates are output to a file
diagnostics = convert_to_bool(obs_par['direct_image_diagnostics'])		#makes diagnostic plot if true

files = glob.glob(os.path.join(path, "*_ima.fits"))				#gets list of filenames in directory

#selects the boundaries of the region where the 2d Gaussian is fit,
#these were selected by eye
rmin = int(obs_par['di_rmin'])
rmax = int(obs_par['di_rmax'])
cmin = int(obs_par['di_cmin'])
cmax = int(obs_par['di_cmax'])

#print "tstart.txt needs to be modified."
#print "It should contain a time before the start time of each VISIT (not orbit), and have 99999. appended to it"

if output ==True: 
	f = open('config/xrefyref.txt', 'w')						#opens file to store positions of reference pixels
	tstart = open('config/tstart.txt', 'w')					#opens file to store exp start times, rounded down


t_last = 0.
for i in enumerate(files):
	ima = pyfits.open(files[i[0]])
	if(ima[0].header['filter'] == obs_par['FILTER']):
		print("filename", files[i[0]])
		LTV1 = ima[1].header['LTV1']					#X offset to get into physical pixels
		LTV2 = ima[1].header['LTV2']					#Y offset to get to physical pixels
		nrow = len(ima[1].data[:,0])
		ncol = len(ima[1].data[0,:])
		t = ima[0].header['expstart']

		dat = ima[1].data[rmin:rmax, cmin:cmax]				#cuts out stamp around the target star
		err = ima[2].data[rmin:rmax, cmin:cmax]

		results = gaussfitter.gaussfit(dat, err)
			
		if diagnostics==True: 
			plt.title("Direct image")
			plt.imshow(dat*ima[0].header['exptime'], origin ='lower',vmin=0, vmax=5000)
			plt.plot(results[2], results[3],marker='x', color='orange', markeredgewidth=3., ms=10, label='centroid', linestyle="none")
			plt.legend(numpoints=1)
			plt.colorbar()
			plt.show()
		
		if output==True: 
			print(t, results[3]+rmin-LTV1, results[2]+cmin-LTV2, file=f)
			if t - t_last > 10./24.:							#checks if this image corresponds to a new visit
				print(np.floor(t), file=tstart)	
				t_last = t
		
		print(t, results[3]+rmin-LTV1, results[2]+cmin-LTV2, files[i[0]]) 		#fit results

if output==True: 
	f.close()
	print(99999., file=tstart)
	tstart.close()
