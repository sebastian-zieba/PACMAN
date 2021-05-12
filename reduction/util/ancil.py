from astropy.io import ascii, fits
import math
import numpy as np

class AncillaryData:
	def __init__(self, obs_par, fit_par=None):
		#used in 0.py and 1.py

		self.path = obs_par['path']

		self.direct_image_output = obs_par['direct_image_output']#flag specifying whether coordinates are output to a file
		self.diagnostics = obs_par['direct_image_diagnostics']#makes diagnostic plot if true

		self.trace_image_output = obs_par['trace_image_output']
		self.all_image_output = obs_par['all_image_output']
		self.trace_image_show = obs_par['trace_image_show']
		self.all_image_show = obs_par['all_image_show']

		# selects the boundaries of the region where the 2d Gaussian is fit,
		# these were selected by eye
		self.rmin = obs_par['di_rmin']
		self.rmax = obs_par['di_rmax']
		self.cmin = obs_par['di_cmin']
		self.cmax = obs_par['di_cmax']
		self.filter = obs_par['FILTER']

		# new ones for 2.py

		#if os.path.isfile("config/filelist.txt"):
		if fit_par != None:

			self.grism = obs_par['GRISM']

			filetable = ascii.read('config/filelist.txt')
			mask = filetable['filter/grism'] == self.grism
			self.files = [self.path + '/' + i for i in filetable['filenames'][mask].data]
			self.t_mjd = filetable['t_mjd'][mask].data
			if 't_bjd' in filetable.keys():
				self.t_bjd = filetable['t_bjd'][mask].data
			self.scans = filetable['scan'][mask].data
			self.orbnum = filetable['norbit'][mask].data
			self.visnum = filetable['nvisit'][mask].data
			self.t_orbit = filetable['t_orbit'][mask].data
			self.t_visit = filetable['t_visit'][mask].data


			refpix = np.genfromtxt("config/xrefyref.txt")  # reads in reference pixels for each visit and sorts them by time
			idx = np.argsort(refpix[:, 0]) #sort by time
			self.refpix = refpix[idx]  # reference pixels from direct image


			self.platescale = 0.13 #IR detector has plate scale of 0.13 arcsec/pixel

			f = fits.open(self.files[0])
			self.POSTARG1 = f[0].header['POSTARG1'] #x-coordinate of the observer requested target offset
			self.POSTARG2 = f[0].header['POSTARG2'] #y-coordinate of the observer requested target offset
			self.LTV1 = int(f[1].header['LTV1'])
			self.LTV2 = int(f[1].header['LTV2'])

			self.BEAMA_i = obs_par['BEAMA_i']  # start of first order trace
			self.BEAMA_f = obs_par['BEAMA_f']  # end of first order trace

			self.subarray_size = f[1].header['SIZAXIS1']  # size of subarray


			self.flat = obs_par['flat']


			self.npix = self.BEAMA_f - self.BEAMA_i  # length of trace



			self.plot_trace = obs_par['plot_trace']
			self.diagnostics = obs_par['diagnostics']
			self.output = obs_par['output']
			if self.output == False: print("NOTE: output is set to False!")

			self.sig_cut = obs_par['sig_cut']
			self.nsmooth = obs_par['nsmooth']
			self.window = obs_par['window']  # window outside of which the background is masked
			self.nvisit = obs_par['nvisit']
			self.norb = obs_par['norb']


			self.skyrmin = obs_par['skyrmin']
			self.skyrmax = obs_par['skyrmax']
			self.skycmin = obs_par['skycmin']
			self.skycmax = obs_par['skycmax']







			self.ra = f[0].header['ra_targ'] * math.pi / 180.0  # stores right ascension
			self.dec = f[0].header['dec_targ'] * math.pi / 180.0  # stores declination

			self.expstart = f[0].header['expstart']  # exposure start time
			self.exptime = f[0].header['exptime']  # exposure time [seconds]

			self.wavegrid = None





			self.one_di_per_visit = obs_par['one_di_per_visit']
			#if self.one_di_per_visit == True: norb = nvisit


			self.torbstart = self.refpix[:, 0]  # start times for each orbit
			self.torbstart = np.append(self.torbstart, 1.0e10)  # appends large value to make get_orbnum routine work

			#self.tstart = tstart

			self.coordtable = []  # table of spacecraft coordinates
			for i in range(max(self.visnum)+1): self.coordtable.append("ancil/bjd_conversion/horizons_results_v" + str(i) + ".txt")


			def make_dict(table):
				return {x['parameter']: x['value'] for x in table}
			fit_par = make_dict(ascii.read(fit_par, Reader=ascii.CommentedHeader))
			self.t0 = fit_par['t0']
			self.period = fit_par['per']
			self.fix_ld = obs_par['fix_ld']
