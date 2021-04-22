import sys
sys.path.insert(0, './models')
import numpy as np
from formatter import FormatParams
from functions import Functions

def calc_astro(t, params, data, funcs, visit):

    flux = np.ones_like(t)
    for i, f in enumerate(funcs.astro): 
        #selects parameters to pass to function
        funcparams = [params[j:j + data.nvisit] for j in funcs.astro_porder[i]]
        flux *= f(t, data, funcparams, visit) 

    return flux 

def calc_sys(t, params, data, funcs, visit):

    flux = np.ones_like(t)
    for i, f in enumerate(funcs.sys): 
        #selects parameters to pass to function
        funcparams = [params[j: j + data.nvisit] for j in funcs.sys_porder[i]]
        flux *= f(t, data, funcparams, visit) 

    return flux 

class Model:
    """
    Stores model fit and related parameters
    """
    def __init__(self, data, myfuncs):
        npoints = len(data.time)

        self.model = np.zeros(npoints)
        self.model_sys = np.zeros(npoints)
        self.model_astro = np.zeros(npoints)
        self.norm_flux = np.zeros(npoints)
        self.phase = np.zeros(npoints)
        self.resid = np.zeros(npoints)
        self.norm_resid = np.zeros(npoints)
        self.chi2 = 0.
        self.chi2red = 0.
        self.rms = 0.
        self.rms_predicted = 1.0e6*np.sqrt(np.mean((data.err/data.flux)**2))
        #self.rms_predicted = 1.0e6*np.sqrt(np.mean(np.sqrt((1./data.flux))**2))                #wong because we binned over several pixels
        self.ln_like = 0.
        self.bic = 0.
        self.params = []
        self.myfuncs = Functions(data, myfuncs)


    def fit(self, data, params):
        #loop over each observation
        for visit in range(data.nvisit):
            #FIXME don't do this every time fit is run
            ind = data.vis_num == visit     

            #t = data.time[ind]
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
        self.ln_like = (-0.5*(np.sum((self.resid/data.err)**2 
            + np.log(2.0*np.pi*(data.err)**2)))
        )
        self.bic = -2.*self.ln_like + data.nfree_param*np.log(data.npoints)
        return self


