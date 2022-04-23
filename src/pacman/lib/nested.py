import numpy as np
import pickle
from datetime import datetime
from scipy.stats import norm
import dynesty
import inspect
from dynesty import plotting as dyplot
import matplotlib.pyplot as plt
import corner
import os
from . import plots
from dynesty import utils as dyfunc


def name_and_args():
    caller = inspect.stack()[1][0]
    args, _, _, values = inspect.getargvalues(caller)
    return [(i, values[i]) for i in args]

def quantile(x, q):                                                             
        return np.percentile(x, [100. * qi for qi in q]) 

def transform_uniform(x,a,b):
    return a + (b-a)*x

def transform_normal(x,mu,sigma):
    return norm.ppf(x,loc=mu,scale=sigma)

def format_params_for_mcmc(params, meta, fit_par):	#FIXME: make sure this works for cases when nvisit>1
    nvisit = int(meta.nvisit)
    theta = []
    meta.fit_par_new = True
    if meta.fit_par_new == False:
        for i in range(len(fit_par)):
                if fit_par['fixed'][i].lower() == "false":
                        if fit_par['tied'][i].lower() == "true": theta.append(params[i*nvisit])
                        else:
                                for j in range(nvisit): theta.append(params[i*nvisit+j])
    else:
        ii = 0
        for i in range(int(len(params)/nvisit)):
                if fit_par['fixed'][ii].lower() == "false":
                        if str(fit_par['tied'][ii]) == "-1":
                            theta.append(params[i*nvisit])
                            ii = ii + 1
                        else:
                            for j in range(nvisit):
                                theta.append(params[i*nvisit+j])
                                ii = ii + 1
                else:
                    ii = ii + 1

    return np.array(theta)



def labels_gen(params, meta, fit_par):
    nvisit = int(meta.nvisit)
    labels = []
    meta.fit_par_new = True
    if meta.fit_par_new == False:
        for i in range(len(fit_par)):
                if fit_par['fixed'][i].lower() == "false":
                        if fit_par['tied'][i].lower() == "true": labels.append(fit_par['parameter'][i])
                        else:
                                for j in range(nvisit): labels.append(fit_par['parameter'][i]+str(j))
    else:
        ii = 0
        for i in range(int(len(params)/nvisit)):
                if fit_par['fixed'][ii].lower() == "false":
                        if str(fit_par['tied'][ii]) == "-1":
                            labels.append(fit_par['parameter'][ii])
                            ii = ii + 1
                        else:
                            for j in range(nvisit):
                                labels.append(fit_par['parameter'][ii]+str(j))
                                ii = ii + 1
                else:
                    ii = ii + 1

    #print('labels', labels)
    return labels





def mcmc_output(samples, params, meta, fit_par, data):	#FIXME: make sure this works for cases when nvisit>1
    nvisit = int(meta.nvisit)
    labels = labels_gen(params, meta, fit_par)

    fig = corner.corner(samples, labels=labels, show_titles=True)
    current_time = datetime.now().time()
    figname = "pairs_"+current_time.isoformat()+".png"
    fig.savefig(figname)


def format_params_for_Model(theta, params, meta, fit_par):
    nvisit = int(meta.nvisit)
    params_updated = []
    iter = 0									#this should be a more informative name FIXME
    meta.fit_par_new = True
    if meta.fit_par_new == False:
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
    else:
        ii = 0
        for i in range(int(len(params)/nvisit)):
                if fit_par['fixed'][ii].lower() == "true":
                        for j in range(nvisit):
                                params_updated.append(params[i*nvisit+j])
                        ii = ii + 1
                else:
                        if str(fit_par['tied'][ii]) == "-1":
                                for j in range(nvisit): params_updated.append(theta[iter])
                                iter += 1
                                ii = ii + 1
                        else:
                                for j in range(nvisit):
                                        params_updated.append(theta[iter])
                                        iter += 1
                                        ii = ii + 1

    return np.array(params_updated)

def nested_sample(data, model, params, file_name, meta, fit_par):
    x = format_params_for_mcmc(params, meta, fit_par)

    ndim = len(x)

    l_args = [params, data, model, meta, fit_par]
    p_args = [data]
    
    #dsampler = dynesty.DynamicNestedSampler(loglike, ptform, ndim,
    #                                        logl_args = l_args,
    #                                       ptform_args = p_args,
    #                                        update_interval=float(ndim))
    #dsampler.run_nested(wt_kwargs={'pfrac': 1.0})#, maxiter = 20000)
    #results = dsampler.results
    sampler = dynesty.NestedSampler(loglike, ptform, ndim, logl_args = l_args,ptform_args = p_args, nlive=meta.run_nlive, bound='single')
    #sampler = dynesty.NestedSampler(loglike, ptform, ndim, logl_args = l_args,ptform_args = p_args, update_interval=float(ndim), nlive=meta.run_nlive, bound='single')
    sampler.run_nested(dlogz=meta.run_dlogz)
    results = sampler.results

    if not os.path.isdir(meta.workdir + meta.fitdir + '/nested_res'):
        os.makedirs(meta.workdir + meta.fitdir + '/nested_res')

    pickle.dump(results, open(meta.workdir + meta.fitdir + '/nested_res/' +  '/nested_out_bin{0}_wvl{1:0.3f}.p'.format(meta.s30_file_counter, meta.wavelength), "wb"))
    results.summary()

    labels = labels_gen(params, meta, fit_par)


    samples, weights = results.samples, np.exp(results.logwt - results.logz[-1])
    mean, cov = dyfunc.mean_and_cov(samples, weights)
    new_samples = dyfunc.resample_equal(samples, weights)

    plots.nested_pairs(new_samples, params, meta, fit_par, data)

    medians = []
    errors_lower = []
    errors_upper = []
    for i in range(ndim):
        q = quantile(new_samples[:, i], [0.16, 0.5, 0.84])
        medians.append(q[1])
        errors_lower.append(abs(q[1] - q[0]))
        errors_upper.append(abs(q[2] - q[1]))

    f_mcmc = open(meta.workdir + meta.fitdir + '/nested_res/' + "/nested_res_bin{0}_wvl{1:0.3f}.txt".format(meta.s30_file_counter, meta.wavelength), 'w')
    for row in zip(errors_lower, medians, errors_upper, labels):
        print('{0: >8}: '.format(row[3]), '{0: >24} '.format(row[1]), '{0: >24} '.format(row[0]), '{0: >24} '.format(row[2]), file=f_mcmc)
    f_mcmc.close()

    updated_params = format_params_for_Model(medians, params, meta, fit_par)
    fit = model.fit(data, updated_params)
    plots.plot_fit_lc2(data, fit, meta, nested=True)
    return medians, errors_lower, errors_upper


def ptform(u, data):
    p = np.zeros_like(u) 
    n = len(data.prior)
    for i in range(n):
        if data.prior[i][0] == 'U':  p[i] = transform_uniform(u[i], 
                                            data.prior[i][1],data.prior[i][2])
        if data.prior[i][0] == 'N':  p[i] = transform_normal(u[i], 
                                            data.prior[i][1],data.prior[i][2])
    return p



def loglike(x, params, data, model, meta, fit_par):
    updated_params = format_params_for_Model(x, params, meta, fit_par)
    fit = model.fit(data, updated_params)
    #print(fit.ln_like)
    return fit.ln_like 
