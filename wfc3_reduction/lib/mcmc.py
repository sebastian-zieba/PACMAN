import emcee
import numpy as np
import pickle
import corner
from datetime import datetime
from scipy.stats import norm
from uncertainties import ufloat
import matplotlib.pyplot as plt
from . import plots


def quantile(x, q):                                                             
        return np.percentile(x, [100. * qi for qi in q]) 

def format_params_for_mcmc(params, meta, fit_par):	#FIXME: make sure this works for cases when nvisit>1
    nvisit = int(meta.nvisit)
    theta = []

    for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "false":
                    if fit_par['tied'][i].lower() == "true": theta.append(params[i*nvisit])
                    else: 
                            for j in range(nvisit): theta.append(params[i*nvisit+j])
    return np.array(theta)


def labels_gen(meta, fit_par):
    nvisit = int(meta.nvisit)
    labels = []

    for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "false":
                    if fit_par['tied'][i].lower() == "true": labels.append(fit_par['parameter'][i])
                    else:
                            for j in range(nvisit): labels.append(fit_par['parameter'][i]+str(j))
    return labels

def mcmc_output(samples, params, meta, fit_par, data):	#FIXME: make sure this works for cases when nvisit>1
    labels = labels_gen(meta, fit_par)

    fig = corner.corner(samples, labels=labels, show_titles=True,quantiles=[0.16, 0.5, 0.84],title_fmt='.4')

    figname = meta.workdir + meta.fitdir + "/pairs_" + meta.fittime + ".png"
    fig.savefig(figname)


def format_params_for_Model(theta, params, meta, fit_par):
    nvisit = int(meta.nvisit)
    params_updated = []
    iter = 0									#this should be a more informative name FIXME
    for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "true": 
                    for j in range(nvisit): 
                            params_updated.append(params[i*nvisit+j])
            else:
                    if fit_par['tied'][i].lower() == "true": 
                            for j in range(nvisit): params_updated.append(theta[iter])
                            iter += 1
                    else: 
                            for j in range(nvisit): 		
                                    params_updated.append(theta[iter])
                                    iter += 1
    return np.array(params_updated)

def mcmc_fit(data, model, params, file_name, meta, fit_par):
    theta = format_params_for_mcmc(params, meta, fit_par)

    ndim, nwalkers = len(theta), meta.run_nwalkers
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (params, data, model, meta, fit_par))

    pos = [theta + 1e-5*np.random.randn(ndim) for i in range(nwalkers)]


    sampler.run_mcmc(pos, meta.run_nsteps, progress=True)
    #sampler.run_mcmc(pos,20000)
    pickle.dump([data, params, sampler.chain], open(meta.workdir + meta.fitdir + "/mcmc_out_"+"{0:0.2f}".format(data.wavelength)+"_"+meta.fittime+".p", "wb"))
    nburn = meta.run_nburn

    if meta.run_nsteps * meta.run_nwalkers > 1000000:
        thin_corner = int((meta.run_nsteps - meta.run_nburn) * meta.run_nwalkers // 100000)
        print('Note: Big Corner plot with many steps. Thinning Plot by factor: {0}'.format(thin_corner))
    else:
        thin_corner = 1

    samples = sampler.chain[:, nburn::thin_corner, :].reshape((-1, ndim))
    mcmc_output(samples, params, meta, fit_par, data)

    samples = sampler.chain[:, nburn:, :].reshape((-1, ndim))
    labels = labels_gen(meta, fit_par)

    plots.plot_chains(ndim, sampler, 0, labels, meta)
    plots.plot_chains(ndim, sampler, nburn, labels, meta)


    medians = []
    errors = []
    errors_new = []
    for i in range(len(theta)):
            q = quantile(samples[:, i], [0.16, 0.5, 0.84])
            medians.append(q[1])
            errors.append(q[2] - q[1])
            errors_new.append(abs((q[2] + q[0])/2 -q[1]))

    f_mcmc = open(meta.workdir + meta.fitdir + "/mcmc_res_{0}.txt".format(meta.fittime), 'w')
    for row in zip(labels, medians, errors_new):
        #print("{: <8} +/-  {: <8}".format(*row))
        xu = ufloat(row[1], row[2])
        print('{0: >12}: '.format(row[0]), '{:13.1uS}'.format(xu), file=f_mcmc)
    f_mcmc.close()

    updated_params = format_params_for_Model(medians, params, meta, fit_par)
    fit = model.fit(data, updated_params)
    print(fit.rms)

    plots.plot_fit_lc(data, fit, meta, mcmc=True)

    return data.wavelength, medians[0], errors[0], samples


def lnprior(theta, data):
    lnprior_prob = 0.
    n = len(data.prior)
    for i in range(n):
        if data.prior[i][0] == 'U': 
            if np.logical_or(theta[i] < data.prior[i][1], 
              theta[i] > data.prior[i][2]): lnprior_prob += - np.inf
        if data.prior[i][0] == 'N': 
            lnprior_prob -= 0.5*(np.sum(((theta[i] - 
              data.prior[i][1])/data.prior[i][2])**2 + 
              np.log(2.0*np.pi*(data.prior[i][2])**2)))
    return lnprior_prob
    


def lnprob(theta, params, data, model, meta, fit_par):
    updated_params = format_params_for_Model(theta, params, meta, fit_par)
    fit = model.fit(data, updated_params)
    lp = lnprior(theta, data)
    return fit.ln_like + lp
