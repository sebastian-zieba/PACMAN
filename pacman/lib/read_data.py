import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import itertools


class Data:
    """
    Reads in and stores raw light curve data
    Args:
        data_file
        obs_par
        fit_par
    """
    def __init__(self, data_file, meta, fit_par, clip_idx=[]):
        def format_prior_for_mcmc(self, meta, fit_par):
                nvisit = int(meta.nvisit)
                prior = []

                for i in range(len(fit_par)):
                    if fit_par['fixed'][i].lower() == "false":
                        prior.append([fit_par['prior'][i], float(fit_par['p1'][i]),
                            float(fit_par['p2'][i])])

                print('prior, read_data', prior)
                return prior

	#read in data and sort by time
        # = np.genfromtxt(data_file)
        d = ascii.read(data_file)
        #d = d[:104]
        #d = d[np.argsort(d[:,5])]


        #removes first exposure from each orbit
        ind = np.diff(d['t_bjd']) < 30./60./24.	# first exposure after 30 mins break
        d = d[1:][ind]

        orb_num = d['iorbit']
        vis_num = d['ivisit']
        t_vis = d['t_visit']
        t_orb = d['t_orbit']

        n = len(d)

        #correct t_obs so that t_obs for the first exposure in orbit = 0
        t_orb_starts = np.zeros(n)
        t_orb_starts[0] = t_orb[0]
        for i, t_orb_diff in enumerate(np.diff(t_orb)):
            if t_orb_diff >= 0: # as long as t_orb increases, use same t_orb_start
                t_orb_starts[i+1] = t_orb_starts[i]
            else:
                t_orb_starts[i+1] = t_orb[i + 1]

        t_orb = t_orb - t_orb_starts
        t_vis = t_vis - min(t_vis)

        #t_delay will =1 if it's the first orbit in a visit. Otherwise =0
        t_delay = np.zeros(n)

        nvisit = int(meta.nvisit)
        norbit = int(meta.norbit)


        
        #orbs_per_visit = norbit/nvisit


        
        #####################################
        #remove first orbit of each visit
        #norbit -= 4
        #TODO: THIS HAS TO BE FIXED!!!
        #TODO: Easiest solution: Do not only save comulative orbit number but also orbit number in visit.
        #TODO: OR EASIER: Save the number of orbits per visit!!
        ind = (orb_num%4 == 0)
        #ind = (orb_num%10 == 0)
        orb_num = orb_num[~ind]
        vis_num = vis_num[~ind]
        t_vis = t_vis[~ind]
        t_vis = t_vis - min(t_vis)
        t_orb = t_orb[~ind]
        t_delay = t_delay[~ind]
        d = d[~ind]

        #TODO: NEEDS FIX!! WONT WORK IF VISIT HASNT 4 ORBITS!!
        ind = (orb_num%4 == 1)
        #ind = (orb_num%10 == 1)
        t_delay[ind] = 1.
        """ind = (orb_num==0)|(orb_num == 6)|(orb_num == 12)|(orb_num == 18)
        t_delay[ind] = 1."""

        #TODO: Will break when user did not use optimal extraction!
        err = np.sqrt(d['var_opt'])
        flux = d['spec_opt']
        time  = d['t_bjd']
        scan_direction = d['scan']

        #FIXME SZ If white it should take wavelength integrated limb darking from file not just at 1.4 micron
        #wavelength = d[0,3]
        if meta.run_fit_white:
            wavelength = 1.4
        else:
            wavelength = d['wave'][0]
        #print "setting wavelength by hand to fix_ld for white lc"


        #fixes limb darkening if "fix_ld" parameter is set to True in obs_par.txt
        if meta.fix_ld == True:
            ld = np.genfromtxt(meta.ld_file)
            
            i = 0
            while(wavelength > ld[i,1]): i += 1
            
            u1 = ld[i, 3]
            u2 = ld[i, 4] 
            fit_par['value'][np.where(fit_par['parameter']=='u1')] = u1
            fit_par['fixed'][np.where(fit_par['parameter']=='u1')] = "true"
            fit_par['value'][np.where(fit_par['parameter']=='u2')] = u2
            fit_par['fixed'][np.where(fit_par['parameter']=='u2')] = "true"

        nfree_param = 0

        for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "false":
                nfree_param += 1

        # The following is for removal of sigma clipped data in the light curve
        # Running this part the first time; clip_idx will always be emtpy
        # After performing the least sq it will decide to clip or not and then possible run this here again
        idx_array = np.arange(len(time), dtype=int)
        #print('readdata:', clip_idx)
        if len(clip_idx) == 0: clip_mask = idx_array
        else:clip_mask = np.bitwise_and.reduce([idx_array!=i for i in clip_idx])

        self.time = time[clip_mask]
        #print(len(time[clip_mask]))
        self.flux = flux[clip_mask]
        print('median log10 raw flux:', np.log10(np.median(flux[clip_mask])))
        self.err = err[clip_mask]
        self.wavelength = wavelength
        self.exp_time = np.median(np.diff(time))*60*60*24 #for supersampling in batman.py
        self.toffset = float(meta.toffset)
        self.nvisit = nvisit
        self.vis_num = vis_num[clip_mask]
        #print(self.vis_num)
        self.orb_num = orb_num[clip_mask]
        self.scan_direction = scan_direction[clip_mask]
        self.t_vis = t_vis[clip_mask]
        self.t_orb = t_orb[clip_mask]
        self.t_delay = t_delay[clip_mask]
        self.parnames = remove_dupl(fit_par['parameter'])
        #print('self.parnames ', self.parnames )
        par_order = {line: i for i, line in enumerate(self.parnames)}
        #print('par_order', par_order)
        self.par_order = par_order
        self.nfree_param = nfree_param
        self.npoints = len(self.time)
        self.dof = self.npoints  - nfree_param
        #TODO: Q: Any need for lc_type?
        self.lc_type = meta.lc_type
        self.all_sys = None
        self.u1 = 0.
        self.u2 = 0.
        print('nfree_param: ', nfree_param)
        #plt.plot(self.t_vis, self.flux, '.k')
        #plt.show()

        self.prior = format_prior_for_mcmc(self, meta, fit_par)

        self.vis_idx = []
        for i in range(nvisit): self.vis_idx.append(self.vis_num == i)
        #print(self.vis_idx)
        #FIXME
        #self.white_systematics = np.genfromtxt("white_systematics.txt")


def remove_dupl(seq):
    #https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
