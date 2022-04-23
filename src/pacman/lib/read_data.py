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

                #print('prior, read_data', prior)
                return prior

	    #read in data
        d = ascii.read(data_file)

        iorbit_sp = meta.iorbit_sp

        #####################################
        # Removes first exposure from each orbit
        # The index where the orbit changes is saved in meta.new_orbit_idx_sp
        if meta.remove_first_exp:
            #removes first exposure from each orbit
            leave_ind = np.ones(len(d), bool)
            leave_ind[meta.new_orbit_idx_sp] = False
            d = d[leave_ind]
            # Remove first exposures from meta.iorbit_sp too
            iorbit_sp = iorbit_sp[leave_ind]
            print('Removed {0} exposures because they were the first exposures in the orbit.'.format(sum(~leave_ind)))
        else:
            print('Leaving the first exposures in every orbit.')

        #####################################
        # Remove first orbit of each visit
        # The NOT cumulative orbit list is saved in meta.iorbit_sp
        if meta.remove_first_orb:
            #removes first exposure from each orbit
            leave_ind = iorbit_sp != 0
            d = d[leave_ind]
            # Remove first orbit from meta.iorbit_sp too
            iorbit_sp = iorbit_sp[leave_ind]
            print('Removed {0} exposures because they were the first orbit in the visit.'.format(sum(~leave_ind)))
        else:
            print('Leaving the first orbit in every visit.')

        n = len(d)

        # t_delay will = 1 if it's the first orbit in a visit. Otherwise = 0
        t_delay = np.zeros(n)
        if meta.remove_first_orb:
            ind = (iorbit_sp == 1)
            t_delay[ind] = 1.
        else:
            ind = (iorbit_sp == 0)
            t_delay[ind] = 1.

        orb_num = d['iorbit']
        vis_num = d['ivisit']
        t_vis = d['t_visit']
        t_orb = d['t_orbit']
        #TODO: Will break when user did not use optimal extraction!
        err = np.sqrt(d['var_opt'])
        flux = d['spec_opt']
        time  = d['t_bjd']
        scan_direction = d['scan']

        # Correct times so that visit and orbit times start at 0 at the beginning of a new visit or orbit
        t_vis = new_time(t_vis)
        t_orb = new_time(t_orb)

        nvisit = int(meta.nvisit)
        norbit = int(meta.norbit)


        #FIXME SZ If white it should take wavelength integrated limb darking from file not just at 1.4 micron
        #wavelength = d[0,3]
        if meta.s30_fit_white:
            meta.wavelength = 1.4
        else:
            meta.wavelength = d['wave'][0]

        #print "setting wavelength by hand to fix_ld for white lc"

        #fixes limb darkening if "fix_ld" parameter is set to True in obs_par.txt
        if meta.fix_ld == True:
            ld = np.genfromtxt(meta.ld_file)
            
            i = 0
            while(meta.wavelength > ld[i,1]): i += 1
            
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
        self.wavelength = meta.wavelength
        self.exp_time = np.median(np.diff(time))*60*60*24 #for supersampling in batman.py
        self.toffset = float(meta.toffset)
        self.nvisit = nvisit ###
        self.vis_num = vis_num[clip_mask] ###
        self.orb_num = orb_num[clip_mask] ###
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
        self.lc_type = meta.lc_type
        self.all_sys = None
        self.u1 = 0.
        self.u2 = 0.
        print('Number of free parameters: ', nfree_param)

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


def new_time(array):
    """
    This functions makes sure the time in a visit (orbit) starts with 0 when a new visit (orbit) starts.
    """
    time_offset = np.zeros(len(array))
    diff_array = np.diff(array)
    current_offset = 0
    for idx, array_i in enumerate(array):
        if idx == 0:
            current_offset = array_i
            time_offset[0] = current_offset
        else:
            if diff_array[idx - 1] > 0:
                time_offset[idx] = current_offset
            else:
                current_offset = array_i
                time_offset[idx] = current_offset
    return array - time_offset
