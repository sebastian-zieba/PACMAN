import numpy as np
from astropy.io import ascii

from .get_ld import get_ld


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
            """
            Reads in the priors from the fit_par file.
            """
            nvisit = int(meta.nvisit)
            prior = []

            for i in range(len(fit_par)):
                # TODO: IMPORTANT! Not only check if its a free parameters but also if it's used!!
                if fit_par['fixed'][i].lower() == "false":
                    prior.append([fit_par['prior'][i], float(fit_par['p1'][i]),
                        float(fit_par['p2'][i])])

            return prior

        #read in data
        d = ascii.read(data_file)

        iorbit_sp = meta.iorbit_sp
        iexp_orb_sp = meta.iexp_orb_sp

        #####################################
        # Removes first exposure from each orbit
        # The index where the orbit changes is saved in meta.new_orbit_idx_sp
        if meta.remove_first_exp:
            # Removes first exposure from each orbit
            leave_ind = np.ones(len(d), bool)
            leave_ind[meta.new_orbit_idx_sp] = False
            d = d[leave_ind]
            # Remove first exposures from meta.iorbit_sp too
            iorbit_sp = iorbit_sp[leave_ind]
            iexp_orb_sp = iexp_orb_sp[leave_ind]
            iexp_orb_sp -= 1 #because we removed iexp_orb = 0 and want the first one in an orbit to be 0 again
            print(f'Removed {sum(~leave_ind)} exposures because they were the first exposures in the orbit.')
        else:
            print('Leaving the first exposures in every orbit.')

        #####################################
        # Remove first orbit of each visit
        # The NOT cumulative orbit list is saved in meta.iorbit_sp
        if len(meta.remove_which_orb) == 1:
            if meta.remove_first_orb:
                # Removes first exposure from each orbit
                leave_ind = iorbit_sp != 0
                d = d[leave_ind]
                # Remove first orbit from meta.iorbit_sp too
                iorbit_sp = iorbit_sp[leave_ind]
                iexp_orb_sp = iexp_orb_sp[leave_ind]
                print(f'Removed {sum(~leave_ind)} exposures because they were the first orbit in the visit.')
        elif len(meta.remove_which_orb) != 1 and meta.remove_first_orb:
            if meta.remove_first_orb:
                masks_orb = []
                for i in range(len(meta.remove_which_orb)):
                    masks_orb.append(iorbit_sp != meta.remove_which_orb[i])
                # Removes chosen orbits from each orbit
                leave_ind = np.bitwise_and(*masks_orb)
                d = d[leave_ind]
                # Remove first orbit from meta.iorbit_sp too
                iorbit_sp = iorbit_sp[leave_ind]
                iexp_orb_sp = iexp_orb_sp[leave_ind]
                print(f'Removed {sum(~leave_ind)} exposures because they were the first orbit in the visit.')
        else:
            print('Leaving the first orbit in every visit.')

        n = len(d)

        # t_delay will = 1 if it's the first orbit in a visit. Otherwise = 0
        t_delay = np.zeros(n)
        if meta.remove_first_orb and meta.remove_which_orb == [0]:
            # if only the first orbit [0] is removed
            ind = (iorbit_sp == 1)
            t_delay[ind] = 1.
        elif meta.remove_first_orb and len(meta.remove_which_orb) != 1:
            # if more than one orbit (eg [0,1]) is removed
            max_orbit_remove = max(meta.remove_which_orb) + 1
            ind = (iorbit_sp == max_orbit_remove)
            t_delay[ind] = 1.
        elif not meta.remove_first_orb:
            ind = (iorbit_sp == 0)
            t_delay[ind] = 1.

        orb_num = d['iorbit'].value
        vis_num = d['ivisit'].value
        t_vis = d['t_visit'].value
        t_orb = d['t_orbit'].value
        err = np.sqrt(d['var_opt'].value)
        flux = d['spec_opt'].value
        time  = d['t_bjd'].value
        scan_direction = d['scan'].value

        # Correct times so that visit and orbit times start at 0 at the beginning of a new visit or orbit
        t_vis = new_time(t_vis)
        t_orb = new_time(t_orb)

        nvisit = int(meta.nvisit)
        norbit = int(meta.norbit)

        # TODO: FIXME SZ If white it should take wavelength integrated limb darking from file not just at 1.4 micron
        if meta.s30_fit_white:
            if meta.grism == 'G102':
                meta.wavelength = 1.0
            elif meta.grism == 'G141':
                meta.wavelength = 1.4
        else:
            meta.wavelength = d['wave'].value[0]

        # fixes limb darkening if "fix_ld" parameter is set to True in obs_par.pcf
        # TODO: Not tested!
        if meta.fix_ld == True:
            ld_path = get_ld(meta)
            ld_file = ascii.read(ld_path)
            
            i = 0
            while(meta.wavelength*1e4 > ld_file['wave_mid'][i]): i += 1
            
            u1 = ld_file['u1'][i]
            u2 = ld_file['u2'][i]
            fit_par['value'][np.where(fit_par['parameter']=='u1')] = u1
            fit_par['fixed'][np.where(fit_par['parameter']=='u1')] = "true"
            fit_par['value'][np.where(fit_par['parameter']=='u2')] = u2
            fit_par['fixed'][np.where(fit_par['parameter']=='u2')] = "true"

        nfree_param = 0
        free_parnames = []
        # TODO: IMPORTANT!
        #  CHECK NOT ONLY THAT FIXED == FALSE BUT ALSO IF THE PARAMETER IS PART OF A MODEL WHICH IS ACTUALLY USED
        for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "false":
                nfree_param += 1
                free_parnames.append(fit_par['parameter'][i].lower())

        # The following is for removal of sigma clipped data in the light curve
        # Running this part the first time; clip_idx will always be emtpy
        # After performing the least sq it will decide to clip or not and then possible run this here again
        idx_array = np.arange(len(time), dtype=int)
        if len(clip_idx) == 0: clip_mask = idx_array
        else: clip_mask = np.bitwise_and.reduce([idx_array!=i for i in clip_idx])

        self.s30_myfuncs = meta.s30_myfuncs
        self.time = time[clip_mask]
        self.flux = flux[clip_mask]
        print('median log10 raw flux of full light curve:', np.log10(np.median(flux[clip_mask])))
        self.err = err[clip_mask]
        self.err_notrescaled = err[clip_mask] # will store the original errorbars and wont be rescaled
        self.wavelength = meta.wavelength
        self.ld_model = meta.ld_model
        self.exp_time = np.median(np.diff(time))*60*60*24 #for supersampling in batman.py
        self.toffset = float(meta.toffset)
        self.nvisit = nvisit
        self.vis_num = vis_num[clip_mask]
        self.orb_num = orb_num[clip_mask]
        self.iexp_orb_sp = iexp_orb_sp[clip_mask]
        self.imax = max(iexp_orb_sp[clip_mask]) + 1 # maximum number of exposures in an orbit (needed for cj model)
        print(f'The highest amount of exposures in an orbit is {self.imax}')
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
        self.free_parnames = free_parnames
        print('Number of free parameters: ', nfree_param)
        print('Names of free parameters: ', free_parnames)
        self.npoints = len(self.time)
        self.dof = self.npoints - self.nfree_param
        #self.lc_type = meta.lc_type # not currently used
        self.all_sys = None
        self.prior = format_prior_for_mcmc(self, meta, fit_par)
        self.vis_idx = []
        for i in range(nvisit): self.vis_idx.append(self.vis_num == i)
        if ('divide_white' in meta.s30_myfuncs) and meta.s30_fit_spec:
            self.white_systematics = np.genfromtxt(meta.white_sys_path)
        self.rescale_uncert = meta.rescale_uncert


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
