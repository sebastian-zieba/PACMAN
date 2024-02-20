import sys

import numpy as np

from ..lib.functions import Functions

sys.path.insert(0, './models')


def calc_astro(t, params, data, funcs, visit):
    flux = np.ones_like(t)
    for i, f in enumerate(funcs.astro):
        # selects parameters to pass to function
        funcparams = [params[j:j + data.nvisit] for j in funcs.astro_porder[i]]
        flux *= f(t, data, funcparams, visit)
    return flux


def calc_sys(t, params, data, funcs, visit):
    flux = np.ones_like(t)
    for i, f in enumerate(funcs.sys):
        # selects parameters to pass to function
        funcparams = [params[j: j + data.nvisit] for j in funcs.sys_porder[i]]
        flux *= f(t, data, funcparams, visit)
    return flux


def calc_gp(idx, params, data, resid, funcs, visit):
    flux = np.ones(int(sum(idx)))
    for i, f in enumerate(funcs.gp):
        # selects parameters to pass to function
        funcparams = [params[j + visit] for j in funcs.gp_porder[i]]
        gp, lnlike = f(idx, data, resid, funcparams)
        flux *= gp
    return [flux, lnlike]


class Model:
    """
    Stores model fit and related parameters
    """
    def __init__(self, data, myfuncs):
        npoints = len(data.time)

        self.model = np.zeros(npoints)
        self.model_sys = np.zeros(npoints)
        self.model_astro = np.zeros(npoints)
        if ('gp_sho' in data.s30_myfuncs) or ('gp_matern32' in data.s30_myfuncs):
            self.model_gp = np.zeros(npoints)
        self.norm_flux = np.zeros(npoints)
        self.phase = np.zeros(npoints)
        self.resid = np.zeros(npoints)
        self.norm_resid = np.zeros(npoints)
        self.chi2 = 0.
        self.chi2red = 0.
        self.rms = 0.
        self.rms_predicted = 1.0e6*np.sqrt(np.mean((data.err/data.flux)**2))
        print(f'The predicted rms is {self.rms_predicted:.2f} ppm')
        # self.rms_predicted = 1.0e6*np.sqrt(np.mean(np.sqrt((1./data.flux))**2))
        # wrong because we binned over several pixels
        self.ln_like = 0.
        self.bic = 0.
        #self.bic_alt = 0.
        self.params = []
        self.myfuncs = Functions(data, myfuncs)

    def fit(self, data, params):
        # loop over each observation
        for visit in range(data.nvisit):
            # FIXME don't do this every time fit is run
            ind = data.vis_num == visit

            # t = data.time[ind]
            t = data.time
            per  = params[data.par_order['per']*data.nvisit + visit]
            t0  = params[data.par_order['t0']*data.nvisit + visit] + data.toffset
            self.phase[ind] = (t[ind] - t0)/per - np.floor((t[ind] - t0)/per)
            self.model_sys[ind] = calc_sys(t[ind], params, data, self.myfuncs, visit)
            self.model_astro[ind] = calc_astro(t[ind], params, data, self.myfuncs, visit)

        self.phase[self.phase > 0.5] -= 1.      #centers phase around 0 for transits

        self.params = params
        self.model = self.model_sys*self.model_astro
        self.data_nosys = data.flux/self.model_sys
        self.norm_flux = data.flux/self.model
        self.all_sys = data.flux/self.model_astro
        self.resid = data.flux - self.model
        self.norm_resid = self.resid/data.flux
        self.chi2 = np.sum((self.resid/data.err)**2)
        self.chi2red = self.chi2/data.dof
        self.rms = 1.0e6*np.sqrt(np.mean((self.resid/data.flux)**2))
        if ('gp_sho' not in data.s30_myfuncs) and ('gp_matern32' not in data.s30_myfuncs):
            self.ln_like = (-0.5*(np.sum((self.resid/data.err)**2
                + np.log(2.0 * np.pi) + 2 * np.log(data.err)))
                # + np.log(2.0*np.pi*(data.err)**2)))
            )
            self.bic = -2. * self.ln_like + data.nfree_param * np.log(data.npoints)
            #self.bic_alt = self.chi2 + data.nfree_param * np.log(data.npoints)
        else:
            self.ln_like = 0.
            for visit in range(data.nvisit):
                # FIXME don't do this every time fit is run
                idx = data.vis_num == visit
                t = data.time[idx]

                gp, gp_ln_like = calc_gp(idx, params, data, self.norm_resid, self.myfuncs, visit)
                self.model_gp[idx] = gp
                self.ln_like += gp_ln_like
            self.params = params
            self.model = self.model_sys * self.model_astro * self.model_gp
            self.data_nosys = data.flux / (self.model_sys * self.model_gp)
            self.norm_flux = data.flux / self.model
            self.all_sys = data.flux / self.model_astro
            self.resid = data.flux - self.model
            self.norm_resid = self.resid / data.flux
            self.chi2 = np.sum((self.resid / data.err) ** 2)
            self.chi2red = self.chi2 / data.dof
            self.rms = 1.0e6 * np.sqrt(np.mean((self.resid / data.flux) ** 2))
            self.bic = -2. * self.ln_like + data.nfree_param * np.log(data.npoints)
            #self.bic_alt = self.chi2 + data.nfree_param * np.log(data.npoints)
        if ('uncmulti' in data.s30_myfuncs) or (data.rescale_uncert): #(data.err[0] != data.err_notrescaled[0]):# or (meta.rescale_uncert): #old:hack so that we dont have to pass meta
            # We want to use err_notrescaled for the chi2_red and BIC calculation
            # because the errors were scaled if one of the following applies:
            # 1) the errorbars were rescaled so that chi2_red = 1
            # 2) uncmulti is turned on and the errorbars were scaled at every step of the sampler
            self.chi2_notrescaled = np.sum((self.resid / data.err_notrescaled) ** 2)
            self.chi2red_notrescaled = self.chi2_notrescaled / (data.dof) # is uncmulti actually a free parameter?
            self.ln_like_notrescaled = (-0.5*(np.sum((self.resid/data.err_notrescaled)**2
                + np.log(2.0 * np.pi) + 2 * np.log(data.err_notrescaled))))
            self.bic_notrescaled = -2. * self.ln_like_notrescaled + data.nfree_param * np.log(data.npoints)
            #self.bic_alt_notrescaled = self.chi2_notrescaled + data.nfree_param * np.log(data.npoints)

        return self
