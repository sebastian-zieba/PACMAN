# This code computes the mean position of the direct image for each visit
import sys
sys.path.insert(0, './util')

from pylab import *
import gaussfitter, ancil, os, glob
from astropy.io import fits
import yaml

yaml_path = "config/obs_par.yaml"
with open(yaml_path, 'r') as file:
	params = yaml.safe_load(file)

ancil = ancil.AncillaryData(params)

files = glob.glob(os.path.join(ancil.path, "*_ima.fits"))				#gets list of filenames in directory


#print "tstart.txt needs to be modified."
#print "It should contain a time before the start time of each VISIT (not orbit), and have 99999. appended to it"

if ancil.output ==True:
	f = open('config/xrefyref.txt', 'w')						#opens file to store positions of reference pixels
	tstart = open('config/tstart.txt', 'w')					#opens file to store exp start times, rounded down

# If file exists, a previous runs has already been executed
# Then, a file called "filelist.txt" should include all _ima.fits file names and if they are a spec or di
prevrun = os.path.exists('config/filelist.txt')

if not prevrun:
	filelist = open('config/filelist.txt', 'w')

	expstart = np.zeros(len(files))
	filtergrism = np.zeros(len(files), dtype=object)

	# Will create a table with the properties of all _ima.fits files at the first run
	for i in enumerate(files):
		filename = i[1]
		ima = fits.open(i[1])
		filtergrism[i[0]] = str(ima[0].header['filter'])
		expstart[i[0]] = ima[0].header['expstart']

	print(expstart)
	print(filtergrism)

	torder = np.argsort(expstart)
	files_sorted = np.array([i.split('/')[-1] for i in files])[torder]
	times = expstart[torder]
	#np.savetxt(filelist, list(zip(files_sorted, filtergrism[torder], expstart[torder])), fmt="%s")

	# Identify orbits and visits

	norbits = []
	nvisits = []
	norbit = 0
	nvisit = 0
	torbs = []
	torb = 0
	norbits.append(norbit)
	nvisits.append(nvisit)
	torbs.append(torb)
	current_first = 0

	for i in enumerate(times[1:]):
		torb = (times[i[0]+1]-times[current_first])*24*60 # time since first exposure in orbit
		print(torb)
		# if two exposures less than 10mins apart -> same orbit
		if np.diff(times)[i[0]]*24*60 < 10:
			torbs.append(torb)
			norbits.append(norbit)
			nvisits.append(nvisit)
		# if two exposures arent in the same orbit and more than an orbital period apart -> not subsequent orbits but a new visit
		elif np.diff(times)[i[0]]*24*60 > 100:
			torb = 0
			torbs.append(torb)
			current_first = i[0]+1
			norbit = 0
			nvisit += 1
			norbits.append(norbit)
			nvisits.append(nvisit)
		# if two exposure more than 10 min apart but less than an orbital period -> subsequent orbits
		else:
			torb = 0
			torbs.append(torb)
			current_first = i[0]+1
			norbit += 1
			norbits.append(norbit)
			nvisits.append(nvisit)

	np.savetxt(filelist, list(zip(files_sorted, filtergrism[torder], nvisits, norbits, torbs, (times-times[0])*24*60)), fmt="%s")
	files = [ancil.path + '/' + i for i in files_sorted[filtergrism[torder] == ancil.FILTER]]
else:
	files, fg,_,_,_,_ = np.loadtxt('config/filelist.txt', dtype=str).T
	files = [ancil.path + '/' + i for i in files[fg == ancil.FILTER]]



t_last = 0.
for i in enumerate(files):
	ima = fits.open(files[i[0]])
	print(ima[0].header['OBSTYPE'])

	if(ima[0].header['filter'] == ancil.FILTER):
		print("filename", files[i[0]])
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
			if t - t_last > 10./24.:							#checks if this image corresponds to a new visit
				print(np.floor(t), file=tstart)	
				t_last = t
		
		print(t, results[3]+ancil.rmin-LTV1, results[2]+ancil.cmin-LTV2, files[i[0]]) 		#fit results
	ima.close()

if ancil.output==True:
	f.close()
	print(99999., file=tstart)
	tstart.close()
