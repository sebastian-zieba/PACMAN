import numpy as np
from . import mpfit
from . import plots
from .plots import plot_raw, plot_fit_lc, plot_fit_lc2
from .formatter import PrintParams
import pickle
from astropy.stats import sigma_clip


def residuals(params, data, model, fjac=None):			
    fit = model.fit(data, params)
    return [0, fit.resid/data.err]

def lsq_fit(fit_par, data, meta, model, myfuncs, noclip=False):
    nvisit = data.nvisit 
    npar = len(fit_par)*nvisit

    #initializes least squares fit parameters
    parinfo = [{'value':0, 'fixed':0, 'limited':[0,0,], 'limits':[0.0,0.0], 
                'step':0.0} for j in range(npar)]
    params_s = []

    #loops through parameters and visits
    #sets initial guess, step size, tie, bounds
    for i in range(int(npar/nvisit)):						
        for j in range(nvisit):						
            parinfo[i*nvisit+j]['value'] = fit_par['value'][i]	
            parinfo[i*nvisit+j]['step'] = fit_par['step_size'][i]
            parinfo[i*nvisit+j]['fixed'] = fit_par['fixed'][i].lower() == "true"
            if j>0 and fit_par['tied'][i].lower() == "true":
                parinfo[i*nvisit+j]['tied'] = 'p[{0}]'.format(nvisit*i)	
            if fit_par['lo_lim'][i].lower() == "true": 		
                parinfo[i*nvisit+j]['limited'][0] = True
                parinfo[i*nvisit+j]['limits'][0] = fit_par['lo_val'][i]
            if fit_par['hi_lim'][i].lower() == "true": 			
                parinfo[i*nvisit+j]['limited'][1] = True
                parinfo[i*nvisit+j]['limits'][1] = fit_par['hi_val'][i]
            params_s.append(fit_par['value'][i])

    params_s = np.array(params_s)

    
    if meta.run_plot_raw_data: plot_raw(data, meta)
    fa = {'data':data, 'model':model}

    # FIXME SZ: NO FILE "white_systematics.txt"
    if meta.run_divide_white:
            sys_vector = np.genfromtxt("white_systematics.txt")
            data.all_sys = sys_vector
            #data.nfree_param -= 2
            #data.dof += 2
#		print "subtracting 2 from dof for divide-white"


    m = mpfit.mpfit(residuals, params_s, functkw=fa, parinfo = parinfo, quiet=True, maxiter=500)

    if noclip == False:
        if sum(np.ma.getmask(sigma_clip(model.resid, sigma=meta.run_clipsigma, maxiters=1))) == 0:
            clip_idx = []
            if meta.run_show_plot: plot_fit_lc2(data, model, meta)
        else:
            clip_idx = np.where(np.ma.getmask(sigma_clip(model.resid, sigma=meta.run_clipsigma, maxiters=1))==True)[0]
            print('Outlier Identified: ', len(clip_idx))
            print('Outlier idx: ', clip_idx)
            if meta.run_show_plot: plot_fit_lc(data, model, meta)

    if noclip == True:
        if meta.run_show_plot: plot_fit_lc2(data, model, meta)

    #print(m.params)

    #rescale error bars based on chi2
    #print("rescaling error bars to get chi2red = 1")
    #print("scale factor = ", np.sqrt(model.chi2red))
    #print data.wavelength, np.sqrt(model.chi2red)
    #data.err = data.err*np.sqrt(model.chi2red)
    #m = mpfit.mpfit(residuals, params_s, functkw=fa, parinfo = parinfo, quiet=True)
    #model = Model(m.params, data, flags)
    
    if m.errmsg: print("MPFIT error message", m.errmsg)

    if meta.run_output:
        f = open(meta.run_out_name, "a")
        print("{0:0.3f}".format(data.wavelength), \
                  "{0:0.6f}".format(m.params[data.par_order['rp']*nvisit]), \
                  "{0:0.6f}".format(m.perror[data.par_order['rp']*nvisit]),\
                  "{0:0.6f}".format(m.params[data.par_order['rp']*nvisit + 1]), \
                  "{0:0.6f}".format(m.perror[data.par_order['rp']*nvisit + 1]),\
                  "{0:0.2f}".format(model.chi2red), file=f)
                   
        #pickle.dump([data, model], open("white_lc_fit.p", "wb"))
        pickle.dump([data, model], open("lsq_fit_" + "{0:0.4f}".format(data.wavelength)+".p", "wb"))
        f.close()

                                                                                                                                                                                                      
    if meta.run_verbose:
        #print "{0:0.3f}".format(data.wavelength), "{0:0.2f}".format(bestfit.chi2red)
        #print data.wavelength, "{0:0.3f}".format(m.params[data.par_order['A1']*nvisit])
        PrintParams(m, data)

    #if meta.run_show_plot: plot_fit_lc(data, model, meta)
    plots.rmsplot(model, data, meta)

    #model = Model(data , myfuncs)

    if noclip == False:
        return data, model, m.params, clip_idx

    if noclip == True:
        return data, model, m.params


