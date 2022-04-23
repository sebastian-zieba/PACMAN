import numpy as np
import pickle
from datetime import datetime
from scipy.stats import norm
import dynesty
import inspect
from dynesty import plotting as dyplot
import os
from . import plots
from dynesty import utils as dyfunc
from . import util


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


def nested_sample(data, model, params, file_name, meta, fit_par):
    theta = util.format_params_for_sampling(params, meta, fit_par)
    #print(theta)
    ndim = len(theta)
    print(ndim)
    l_args = [params, data, model, meta, fit_par]
    p_args = [data]
    
    #dsampler = dynesty.DynamicNestedSampler(loglike, ptform, ndim,
    #                                        logl_args = l_args,
    #                                       ptform_args = p_args,
    #                                        update_interval=float(ndim))
    #dsampler.run_nested(wt_kwargs={'pfrac': 1.0})#, maxiter = 20000)
    #results = dsampler.results

    print('Run dynesty...')
    sampler = dynesty.NestedSampler(loglike, ptform, ndim, logl_args = l_args,ptform_args = p_args,update_interval=float(ndim), nlive=meta.run_nlive)
    sampler.run_nested(dlogz=meta.run_dlogz)
    results = sampler.results

    if not os.path.isdir(meta.workdir + meta.fitdir + '/nested_res'):
        os.makedirs(meta.workdir + meta.fitdir + '/nested_res')

    pickle.dump(results, open(meta.workdir + meta.fitdir + '/nested_res/' +  '/nested_out_bin{0}_wvl{1:0.3f}.p'.format(meta.s30_file_counter, meta.wavelength), "wb"))
    results.summary()

    labels = meta.labels
    #labels = labels_gen(params, meta, fit_par)


    samples, weights = results.samples, np.exp(results.logwt - results.logz[-1])
    mean, cov = dyfunc.mean_and_cov(samples, weights)
    new_samples = dyfunc.resample_equal(samples, weights)


    # Plot a summary of the run.
    #rfig, raxes = dyplot.runplot(results)
    #plt.savefig(meta.workdir + meta.fitdir + '/nested_runplot_' + meta.fittime + '.png')
    # Plot traces and 1-D marginalized posteriors.
    #tfig, taxes = dyplot.traceplot(results)
    #plt.savefig(meta.workdir + meta.fitdir + '/nested_traceplot_' + meta.fittime + '.png')
    # Plot the 2-D marginalized posteriors.
    #cfig, caxes = dyplot.cornerplot(results, show_titles=True, title_fmt='.4',labels=labels, color='blue', hist_kwargs=dict(facecolor='blue', edgecolor='blue'))
    #plt.savefig(meta.workdir + meta.fitdir + '/nested_res/' +  '/nested_cornerplot_' + meta.fittime + '.png')
    #plt.show()

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

    updated_params = util.format_params_for_Model(medians, params, meta, fit_par)
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
    updated_params = util.format_params_for_Model(x, params, meta, fit_par)
    fit = model.fit(data, updated_params)
    return fit.ln_like 
