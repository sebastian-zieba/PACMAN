import numpy as np
import matplotlib.pyplot as plt

class Data:
    """
    Reads in and stores raw light curve data
    Args:
        data_file
        obs_par
        fit_par
    """
    def __init__(self, data_file, obs_par, fit_par):
        def format_prior_for_mcmc(self, obs_par, fit_par):
                nvisit = int(obs_par['nvisit'])				
                prior = []

                for i in range(len(fit_par)):
                    if fit_par['fixed'][i].lower() == "false":
                        if fit_par['tied'][i].lower() == "true": 
                            prior.append([fit_par['prior'][i], float(fit_par['p1'][i]), 
                                float(fit_par['p2'][i])])
                        else: 
                            for j in range(nvisit): 
                                prior.append([fit_par['prior'][i], float(fit_par['p1'][i]), 
                                    float(fit_par['p2'][i])])
                return prior

	#read in data and sort by time
        d = np.genfromtxt(data_file)
        d = d[np.argsort(d[:,5])]   #FIXME (put indices in a file, or add header)


        #removes first exposure from each orbit
        ind = np.diff(d[:,5]) < 30./60./24.	# first exposure after 30 mins break
        d = d[1:][ind]

        orb_num = d[:,7]
        vis_num = d[:,6]
        t_vis = d[:,9]
        t_orb = d[:,10]


        n = len(d)




        t_orb_starts = np.zeros(n)
        ind = np.diff(d[:,5]) < 30./60./24.
        t_orb_startsi = np.concatenate(([t_orb[0]], t_orb[1:][~ind]))

        import itertools
        Y = [(x, len(list(y))) for x, y in itertools.groupby(orb_num)]
        #print(Y)

        t_orb_starts = (np.array([[t_orb_startsi[i]]*Y[i][1] for i in range(len(Y))])).flatten()

        t_orb = t_orb - t_orb_starts
        t_vis = t_vis - min(t_vis)

        t_delay = np.zeros(n)

        nvisit = int(obs_par['nvisit'])
        norbit = int(obs_par['norb'])


        
        #orbs_per_visit = norbit/nvisit


        
        #####################################
        #remove first orbit of each visit
        #norbit -= 4
        ind = (orb_num == 0)

        orb_num = orb_num[~ind]
        vis_num = vis_num[~ind]
        t_vis = t_vis[~ind]
        t_vis = t_vis - min(t_vis)
        t_orb = t_orb[~ind]
        t_delay = t_delay[~ind]
        d = d[~ind]

        ind = (orb_num==1)
        t_delay[ind] = 1.
        """ind = (orb_num==0)|(orb_num == 6)|(orb_num == 12)|(orb_num == 18)
        t_delay[ind] = 1."""

        err = np.sqrt(d[:,2])
        flux = d[:,1]
        time  = d[:,5]
        scan_direction = d[:,8]

        wavelength = d[0,3]
        #wavelength = 1.4
        #print "setting wavelength by hand to fix_ld for white lc"


        #fixes limb darkening if "fix_ld" parameter is set to True in obs_par.txt
        if obs_par['fix_ld'] == True:
            ld = np.genfromtxt(obs_par['ld_file'])
            
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
                if fit_par['tied'][i].lower() == "true": nfree_param += 1
                else: nfree_param += nvisit

        self.time = time
        self.flux = flux
        self.err = err
        self.wavelength = wavelength
        self.exp_time = np.median(np.diff(time))*60*60*24 #float(obs_par['exp_time'])
        self.toffset = float(obs_par['toffset'])
        self.nvisit = nvisit
        self.vis_num = vis_num
        self.orb_num = orb_num
        self.scan_direction = scan_direction
        self.t_vis = t_vis
        self.t_orb = t_orb
        self.t_delay = t_delay
        self.parnames = fit_par['parameter']
        par_order = {line['parameter']: i for i, line in enumerate(fit_par)}
        self.par_order = par_order
        self.nfree_param = nfree_param
        self.npoints = len(self.time)
        self.dof = self.npoints  - nfree_param
        self.lc_type = obs_par['lc_type']
        self.all_sys = None
        self.u1 = 0.
        self.u2 = 0.

        #plt.plot(self.t_vis, self.flux, '.k')
        #plt.show()

        self.prior = format_prior_for_mcmc(self, obs_par, fit_par)

        self.vis_idx = []
        for i in range(nvisit): self.vis_idx.append(self.vis_num == i)

        #FIXME
        #self.white_systematics = np.genfromtxt("white_systematics.txt")