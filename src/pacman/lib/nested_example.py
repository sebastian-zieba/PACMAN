import numpy as np
import matplotlib.pyplot as plt
import spiderman as sp
import batman
from multiprocessing import Pool
from IPython.display import display, Math
import corner
import time
import emcee
import matplotlib as mpl
from scipy.optimize import minimize
from scipy.special import ndtri
from dynesty import NestedSampler
import dynesty
from dynesty import plotting as dyplot
from dynesty.utils import resample_equal
from matplotlib.offsetbox import AnchoredText



magK = 7.55
efficieny = 0.66
texp = 138.381

noise_ref = 103e-6
magK_ref = 10.767
efficieny_ref = 1
texp_ref = 103


def T14(per, ars, rprs, inc_rad):
    b = ars * np.cos(inc_rad)
    return per/np.pi * np.arcsin(1/ars * 1/np.sin(inc_rad) * np.sqrt((1+rprs)**2-b**2))


t0=0.
per=0.224
rprs=0.02046
ars=2.063
inc_rad=90*np.pi/180


D=T14(per, ars, rprs, inc_rad)


def toymodel(t, albedo, redist):
    
    t0=0.
    per=0.224
    rprs=0.02046
    ars=2.063
    inc_rad=90*np.pi/180

    Ts = 4250
    Ms = 0.638 
    a_abs = ((per/365.25)**2*Ms)**(1/3)
    
    
    params = batman.TransitParams()       #object to store transit parameters
    params.t0 = t0                       #time of inferior conjunction
    params.per = per                      #orbital period
    params.rp = rprs                      #planet radius (in units of stellar radii)
    params.a = ars                     #semi-major axis (in units of stellar radii)
    params.inc = inc_rad*180/  np.pi                    #orbital inclination (in degrees)
    params.ecc = 0.                       #eccentricity
    params.w = 90.                        #longitude of periastron (in degrees)
    params.limb_dark = "linear"        #limb darkening model
    params.u = [0.04]      #limb darkening coefficients [u1, u2, u3, u4]

    mt = batman.TransitModel(params, t)    #initializes model
    fluxt = mt.light_curve(params)                    #calculates light curve

    spider_params = sp.ModelParams(brightness_model="Louden")

    spider_params.n_layers= 4
    spider_params.grid_size= 3  
    spider_params.t0= t0             # Central time of PRIMARY transit [days]
    spider_params.per= per       # Period [days]
    spider_params.a_abs= a_abs        # The absolute value of the semi-major axis [AU]
    spider_params.inc= inc_rad*180/  np.pi           # Inclination [degrees]
    spider_params.ecc= 0.0              # Eccentricity
    spider_params.w= 90                 # Argument of periastron
    spider_params.rp= rprs            # Planet to star radius ratio
    spider_params.a= ars             # Semi-major axis scaled by stellar radius
    spider_params.p_u1= 0               # Planetary limb darkening parameter
    spider_params.p_u2= 0               # Planetary limb darkening parameter


    spider_params.T_s = Ts    # Temperature of the star
    spider_params.albedo = albedo    # Temperature of the star
    spider_params.insol= 5.67e-8 * Ts**4/ars**2   # Temperature of the star    
    spider_params.redist= redist
    spider_params.l1 = 1.1e-6       # The starting wavelength in meters
    spider_params.l2 = 1.7e-6       # The ending wavelength in meters
    spider_params.T_int = 0

    spider_params.thermal = True
    
    fluxsp = spider_params.lightcurve(t)
    
    return fluxt*fluxsp


t_model_offset = 0.05
t_obs_offset = 0.04

t_model = np.linspace(-per/2-D/2-t_model_offset, per/2+D/2+t_model_offset, 1000)
f_model = toymodel(t_model, 0,0)


t_obs = np.linspace(-per/2-D/2-t_obs_offset, per/2+D/2+t_obs_offset, 100)
f_obs = toymodel(t_obs, 0,0)

np.random.seed(69420)

noise = noise_ref *(10**(-0.4*(magK_ref-magK)))**(1/2)*(efficieny_ref/efficieny)**(1/2)
f_noise = np.random.normal(0, noise, len(t_obs))

datat = np.array(t_obs)
dataf = np.array((f_obs+f_noise))
dataf_err = np.array([noise]*len(datat))


fig, ax = plt.subplots(1,1, figsize=(12,5))

ax.plot(t_model, f_model, c='r', lw=3, alpha=0.6, zorder=-3)
ax.errorbar(datat, dataf, yerr=dataf_err, fmt='.')
ax.set_xlabel('time (days)')
ax.set_ylabel('flux')
plt.tight_layout()
#ax.set_xlim(-0.225, 0.67)
#plt.title('observations {0} minutes relative to eclipse'.format(int(delay*1440)))
#plt.savefig('observations{0}.png'.format(int(delay*1440)), dpi=250)



x,y,yerr = datat, dataf, dataf_err





def prior_transform(theta):
    """
    A function defining the tranform between the parameterisation in the unit hypercube
    to the true parameters.

    Args:
        theta (tuple): a tuple containing the parameters.
        
    Returns:
        tuple: a new tuple or array with the transformed parameters.
    """

    aprime, bprime = theta # unpack the parameters (in their unit hypercube form)
    
    amin = 0  # lower bound on uniform prior on c
    amax = 1 #  er bound on uniform prior on c

    bmin = 0  # lower bound on uniform prior on c
    bmax = 0.5   # upper bound on uniform prior on c

    a = aprime*(amax-amin) + amin  # convert back to c
    b = bprime*(bmax-bmin) + bmin  # convert back to c

    return (a, b)


# set the natural logarithm of 2pi, so that it doesn't have to be recalculated
LN2PI = np.log(2.*np.pi)
LNSIGMA = np.log(noise) # natural log of the data noise standard deviation

M = len(x)

def loglikelihood_dynesty(theta):
    """
    The log-likelihood function.
    """

    albedo,redist = theta # unpack the parameters

    # normalisation
    norm = -0.5*M*LN2PI - M*LNSIGMA

    # chi-squared (data, sigma and x are global variables defined early on in this notebook)
    #chisq = np.sum(((dataf-(phasecurve(datat, per, ars, rprs, inc_rad, fpfs, ampl, D)))/(noise))**2)
    chisq = np.sum(((y-(toymodel(x, albedo,redist)))/(noise))**2)
    
    return norm - 0.5*chisq






nlive = 250      # number of live points
bound = 'multi'   # use MutliNest algorithm for bounds
ndims = 2         # two parameters
sample = 'unif'   # uniform sampling
tol = 0.1         # the stopping criterion



sampler = NestedSampler(loglikelihood_dynesty, prior_transform, ndims,
                        bound=bound, sample=sample, nlive=nlive)

t00 = time.time()
sampler.run_nested(dlogz=tol, print_progress=True) # don't output progress bar
t11 = time.time()

timedynesty = (t11-t00)

print("Time taken to run 'dynesty' (in static mode) is {} seconds".format(timedynesty))



res = sampler.results # get results dictionary from sampler

logZdynesty = res.logz[-1]        # value of logZ
logZerrdynesty = res.logzerr[-1]  # estimate of the statistcal uncertainty on logZ

print("log(Z) = {} ± {}".format(logZdynesty, logZerrdynesty))



print(res.summary())

truths = [0,0]

# plot the resulting posteriors
#mpl.rcParams.update({'font.size': 16})

#def plotposts(samples, truths=[m,c]):
def plotposts(samples, **kwargs):
    """
    Function to plot posteriors using corner.py and scipy's gaussian KDE function.
    """
    if "truths" not in kwargs:
        kwargs["truths"] = truths

    fig = corner.corner(samples, labels=[r'$albedo$', r'$redist$'], hist_kwargs={'density': True}, **kwargs,show_titles=True)

    #plotb KDE smoothed version of distributions
    #for axidx, samps in zip([0, 3], samples.T):
    #    kde = gaussian_kde(samps)
    #    xvals = fig.axes[axidx].get_xlim()
    #    xvals = np.linspace(xvals[0], xvals[1], 100)
     #   fig.axes[axidx].plot(xvals, kde(xvals), color='firebrick')
    #fig.savefig('test_dynesty_TOI-2431b.jpeg', dpi=250)
    
    
    
# get function that resamples from the nested samples to give sampler with equal weight


# draw posterior samples
weights = np.exp(res['logwt'] - res['logz'][-1])
samples_dynesty = resample_equal(res.samples, weights)


print('Number of posterior samples is {}'.format(len(samples_dynesty)))

# plot using corner.py
plotposts(samples_dynesty)




# plot initial run (left)
fg, ax = dyplot.cornerpoints(res, cmap='plasma', truths=truths,
                             kde=False)


# plotting the original run
fig, axes = dyplot.traceplot(res, truths=truths, truth_color='black',
                             show_titles=True, title_kwargs={'fontsize': 28, 'y': 1.05},
                             trace_cmap='plasma', kde=False,
                             connect=True, connect_highlight=range(2),
                             fig=plt.subplots(2, 2, figsize=(14, 12)))
fig.tight_layout()


  # initialize figure
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# plot 6 snapshots over the course of the run
for i, a in enumerate(axes.flatten()):
    it = int((i+1)*res.niter/8.)
    # overplot the result onto each subplot
    temp = dyplot.boundplot(res, dims=(0, 1), it=it, prior_transform=prior_transform, max_n_ticks=3,
                            show_live=True, span=[(0,1), (0,0.5)], fig=(fig, a))
    a.set_title('Iteration {0}'.format(it), fontsize=26)
    a.set_xlabel('albedo')
    a.set_ylabel('redist')


anchored_text = AnchoredText(r"niter {0}".format(res.niter), loc=2)
axes[1,2].add_artist(anchored_text)

fig.tight_layout()
#fig.savefig('snaps-AB001-F001-lp300-multi-rand69.jpeg', dpi=250)






