import numpy as np
import pickle
from scipy.stats import norm
import dynesty
import inspect
import os
from . import plots
from dynesty import utils as dyfunc
from . import util


def transform_uniform(x,a,b):
    return a + (b-a)*x


def transform_normal(x,mu,sigma):
    return norm.ppf(x,loc=mu,scale=sigma)


def nested_sample(data, model, params, file_name, meta, fit_par):
    theta = util.format_params_for_sampling(params, meta, fit_par)
    ndim = len(theta)
    l_args = [params, data, model, meta, fit_par]
    p_args = [data]

    print('Run dynesty...')
    if meta.run_dynamic:
        sampler = dynesty.DynamicNestedSampler(loglike, ptform, ndim, logl_args = l_args, ptform_args = p_args,
                                               update_interval=float(ndim), bound=meta.run_bound,
                                               sample=meta.run_sample)
        sampler.run_nested(wt_kwargs={'pfrac': 1.0}, print_progress=True)#, maxiter = 20000)
    else:
        sampler = dynesty.NestedSampler(loglike, ptform, ndim, logl_args = l_args, ptform_args = p_args,
                                        update_interval=float(ndim), nlive=meta.run_nlive, bound=meta.run_bound,
                                        sample=meta.run_sample)
        sampler.run_nested(dlogz=meta.run_dlogz, print_progress=True)

    results = sampler.results

    if not os.path.isdir(meta.workdir + meta.fitdir + '/nested_res'):
        os.makedirs(meta.workdir + meta.fitdir + '/nested_res')

    pickle.dump(results, open(meta.workdir + meta.fitdir + '/nested_res/' +  '/nested_out_bin{0}_wvl{1:0.3f}.p'.format(meta.s30_file_counter, meta.wavelength), "wb"))
    results.summary()

    labels = meta.labels

    samples, weights = results.samples, np.exp(results.logwt - results.logz[-1])
    mean, cov = dyfunc.mean_and_cov(samples, weights)
    new_samples = dyfunc.resample_equal(samples, weights)

    plots.dyplot_runplot(results, meta)
    plots.dyplot_traceplot(results, meta)
    plots.dyplot_cornerplot(results, meta)
    plots.nested_pairs(new_samples, params, meta, fit_par, data)

    medians = []
    errors_lower = []
    errors_upper = []
    for i in range(ndim):
        q = util.quantile(new_samples[:, i], [0.16, 0.5, 0.84])
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
