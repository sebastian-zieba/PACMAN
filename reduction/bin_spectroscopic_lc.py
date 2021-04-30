#This code reads in the optimally extracted lightcurve and bins it into channels 5 pixels wide, following Berta '12
from numpy import *
from pylab import *
from astropy.io import ascii
from scipy import signal
import matplotlib.pyplot as plt
import os

def make_dict(table): return {x['parameter']: x['value'] for x in table}

def weighted_mean(data, err):				#calculates the weighted mean for data points data with std err
	weights = 1.0/err**2.
	mu = np.sum(data*weights)/np.sum(weights)
	var = 1.0/np.sum(weights)
	return [mu, np.sqrt(var)]				#returns weighted mean and variance

#what bins do you want?
#wave_bins = np.linspace(1.125, 1.65, 22)*1e4
wave_bins = np.linspace(1.125, 1.65, 12)*1e4
print(wave_bins)
#wave_bins = np.linspace(1.139, 1.631, 12)*1e4
#wave_bins = np.array([11400, 11800, 12200,12600,13000, 13400, 13800, 14200, 14600, 15000, 15400, 15800, 16200])

#Tom Evan wave bins
#wave_bins = np.array([1.129595435219961752e+00, 1.180470435219961756e+00, 1.231345435219961759e+00, 1.282220435219961763e+00, 1.333095435219961766e+00, 1.383970435219961770e+00, 1.434845435219961773e+00, 1.485720435219961777e+00, 1.536595435219961781e+00, 1.587470435219961784e+00, 1.638345435219961788e+00])*1.e4

#reads in spectra

lst = os.listdir('extracted_lc')
print(lst.sort())
d = np.genfromtxt(lst.sort()[-1] + "/lc_spec.txt")

obs_par = make_dict(ascii.read("config/obs_par.txt", Reader=ascii.CommentedHeader))
nexp = 104#int(obs_par['nexp'])			#number of exposures
npix = 181                                      #width of spectrum in pixels (BEAMA_f - BEAMA_i) 
d = d.reshape(nexp, npix, -1)			#reshapes array by exposure

w = d[0,:, 4]
f = d[0, :, 2]

w_hires = np.linspace(w.min(), w.max(), 10000)
oversample_factor = len(w_hires)/len(w)*1.0

#stores the indices corresponding to the wavelength range in each bin
wave_inds = []
lo_res_wave_inds = []
for i in range(len(wave_bins)- 1): wave_inds.append((w_hires >= wave_bins[i])&(w_hires <= wave_bins[i+1]))
for i in range(len(wave_bins)- 1): lo_res_wave_inds.append((w >= wave_bins[i])&(w <= wave_bins[i+1]))


dirname = "extracted_sp/" + datetime.strftime(datetime.now(), '%Y-%m-%d_%H:%M')
if not os.path.exists(dirname): os.makedirs(dirname)

for i in range(len(wave_bins) - 1):
    wave = (wave_bins[i] + wave_bins[i+1])/2./1.e4
    outname = dirname + "/speclc" + "{0:.3f}".format(wave)+".txt"
    #outname = "wasp33b_" + "{0:.4f}".format(wave)+".txt"
    outfile = open(outname, 'w')
    

    print("#time", '\t\t', "photoelectrons", '\t', "error", '\t\t', "visit", '\t', "orbit", '\t', "scan", '\t', "wave_center", '\t', "wave_start",  '\t', "wave_end", file=outfile)
    for j in range(nexp):
        time, phase, visnum, orbnum, scan = d[j,0,0], d[j,0,1], d[j,0,5], d[j,0,6], d[j,0,7]

        f_interp = np.interp(w_hires, w, d[j,:,2])
        variance_interp = np.interp(w_hires, w, d[j,:,3])

        
        #accounts for decrease in precision when spectrum is oversampled
        variance_interp *= oversample_factor


        fluxes = f_interp[wave_inds[i]]
        errs = np.sqrt(variance_interp[wave_inds[i]])
    
        meanflux, meanerr = weighted_mean(fluxes, errs)		

        print(phase, meanflux, meanerr**2, wave, 0., time, visnum, orbnum, scan, file=outfile)
            
        #print wave, np.sum(d[j, lo_res_wave_inds[i],2]) 

#print wave, 1.0*sum(wave_inds)/len(w_hires), meanflux, meanerr

plt.plot(w_hires, f_interp)
for wave in wave_bins: plt.axvline(wave, color = '0.5')
plt.ylabel("Photelectrons")
plt.xlabel("Wavelength (angstroms)")
plt.show()
