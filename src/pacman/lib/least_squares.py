import numpy as np
import os
from . import mpfit
from . import plots
from .formatter import PrintParams
import pickle
from astropy.stats import sigma_clip


def residuals(params, data, model, fjac=None):			
    fit = model.fit(data, params)
    return [0, fit.resid/data.err]


def lsq_fit(fit_par, data, meta, model, myfuncs, noclip=False):
    #TODO: noclip = True should be standard
    nvisit = data.nvisit 
    npar = len(data.parnames)*nvisit

    #initializes least squares fit parameters
    parinfo = [{'value':0, 'fixed':0, 'limited':[0,0,], 'limits':[0.0,0.0], 
                'step':0.0} for j in range(npar)]
    params_s = []

    #loops through parameters and visits
    #sets initial guess, step size, tie, bounds

    ii=0
    for i in range(int(len(data.parnames))):
        if str(fit_par['tied'][ii]) == "-1":
            for j in range(nvisit):
                parinfo[i*nvisit+j]['value'] = fit_par['value'][ii]
                parinfo[i*nvisit+j]['step'] = fit_par['step_size'][ii]
                parinfo[i*nvisit+j]['fixed'] = fit_par['fixed'][ii].lower() == "true"
                if j>0 and str(fit_par['tied'][ii]) == "-1":
                    parinfo[i*nvisit+j]['tied'] = 'p[{0}]'.format(nvisit*i)
                if fit_par['lo_lim'][ii].lower() == "true":
                    parinfo[i*nvisit+j]['limited'][0] = True
                    parinfo[i*nvisit+j]['limits'][0] = fit_par['lo_val'][ii]
                if fit_par['hi_lim'][ii].lower() == "true":
                    parinfo[i*nvisit+j]['limited'][1] = True
                    parinfo[i*nvisit+j]['limits'][1] = fit_par['hi_val'][ii]
                params_s.append(fit_par['value'][ii])
            ii = ii+1
        else:
            for j in range(nvisit):
                parinfo[i * nvisit + j]['value'] = fit_par['value'][ii]
                parinfo[i * nvisit + j]['step'] = fit_par['step_size'][ii]
                parinfo[i * nvisit + j]['fixed'] = fit_par['fixed'][ii].lower() == "true"
                if j > 0 and str(fit_par['tied'][ii]) == "-1":
                    parinfo[i * nvisit + j]['tied'] = 'p[{0}]'.format(nvisit * i)
                if fit_par['lo_lim'][ii].lower() == "true":
                    parinfo[i * nvisit + j]['limited'][0] = True
                    parinfo[i * nvisit + j]['limits'][0] = fit_par['lo_val'][ii]
                if fit_par['hi_lim'][ii].lower() == "true":
                    parinfo[i * nvisit + j]['limited'][1] = True
                    parinfo[i * nvisit + j]['limits'][1] = fit_par['hi_val'][ii]
                params_s.append(fit_par['value'][ii])
                ii = ii + 1


    params_s = np.array(params_s)
    if meta.save_raw_lc_plot: plots.plot_raw(data, meta)
    fa = {'data':data, 'model':model}

    if meta.run_divide_white and meta.s30_fit_spec:
        sys_vector = np.genfromtxt(meta.workdir + "/white_systematics.txt")
        data.all_sys = sys_vector
        print("subtracting 2 from dof for divide-white")
        #data.nfree_param -= 2
        #data.dof += 2

    print('\nRuns MPFIT... ')
    m = mpfit.mpfit(residuals, params_s, functkw=fa, parinfo = parinfo, quiet=True)

    if noclip == False:
        #if user wants to sigma clip but there is nothing to clip:
        if sum(np.ma.getmask(sigma_clip(model.resid, sigma=meta.run_clipsigma, maxiters=1))) == 0:
            clip_idx = []
            if meta.save_fit_lc_plot: plots.plot_fit_lc2(data, model, meta)
        else:
            #if user wants to sigma clip and there are outliers:
            clip_idx = np.where(np.ma.getmask(sigma_clip(model.resid, sigma=meta.run_clipsigma, maxiters=1))==True)[0]
            print('Outlier Identified: ', len(clip_idx))
            print('Outlier idx: ', clip_idx)
            if meta.save_fit_lc_plot: plots.plot_fit_lc(data, model, meta)


    #print(m.params)

    #rescale error bars based on chi2
    #print("rescaling error bars to get chi2red = 1")
    #print("scale factor = ", np.sqrt(model.chi2red))
    #print data.wavelength, np.sqrt(model.chi2red)
    #data.err = data.err*np.sqrt(model.chi2red)
    #m = mpfit.mpfit(residuals, params_s, functkw=fa, parinfo = parinfo, quiet=True)
    #model = Model(m.params, data, flags)
    
    if m.errmsg: print("MPFIT error message", m.errmsg)


    # if meta.run_output:
    #     f = open(meta.run_out_name, "a")
    #     print("{0:0.3f}".format(data.wavelength), \
    #               "{0:0.6f}".format(m.params[data.par_order['rp']*nvisit]), \
    #               "{0:0.6f}".format(m.perror[data.par_order['rp']*nvisit]),\
    #               "{0:0.6f}".format(m.params[data.par_order['rp']*nvisit + 1]), \
    #               "{0:0.6f}".format(m.perror[data.par_order['rp']*nvisit + 1]),\
    #               "{0:0.2f}".format(model.chi2red), file=f)
    #     #pickle.dump([data, model], open("white_lc_fit.p", "wb"))
    #     pickle.dump([data, model], open("lsq_fit_" + "{0:0.4f}".format(data.wavelength)+".p", "wb"))
    #     f.close()

    if meta.run_verbose:
        if not os.path.isdir(meta.workdir + meta.fitdir + '/lsq_res'):
            os.makedirs(meta.workdir + meta.fitdir + '/lsq_res')
        f_lsq = open(meta.workdir + meta.fitdir + '/lsq_res/' + "/lsq_res_bin{0}_wvl{1:0.3f}.txt".format(meta.s30_file_counter, meta.wavelength), 'w')
        PrintParams(m, data, savefile=f_lsq)
        PrintParams(m, data)

    if meta.save_allan_plot:
        plots.rmsplot(model, data, meta)

    if noclip == False:
        return data, model, m.params, clip_idx, m

    if noclip == True:
        return data, model, m.params, m
