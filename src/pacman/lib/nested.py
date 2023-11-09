import pickle
from multiprocessing import Pool

import dynesty
import numpy as np
from dynesty import utils as dyfunc
from scipy.stats import norm

from . import plots
from . import util


def transform_uniform(x, a, b):
    """
    Prior transform for uniform priors.
    """
    return a + (b-a)*x


def transform_normal(x, mu, sigma):
    """
    Prior transform for normal priors.
    """
    return norm.ppf(x, loc=mu, scale=sigma)


def nested_sample(data, model, params, file_name, meta, fit_par):
    """
    Calls the dynesty package and does the sampling.
    """
    nvisit = int(meta.nvisit)

    # Create the nested_res directory
    util.create_res_dir(meta)

    # Setting up parameters for sampler
    theta = util.format_params_for_sampling(params, meta, fit_par)
    ndim = len(theta)

    fixed_array = np.array(fit_par['fixed'])
    tied_array = np.array(fit_par['tied'])
    free_array = util.return_free_array(nvisit, fixed_array, tied_array)
    l_args = [params, data, model, nvisit, fixed_array, tied_array, free_array]
    p_args = [data]

    # Setting up multiprocessing
    if hasattr(meta, 'ncpu') and meta.ncpu > 1:
        print('Using multiprocessing...')
        pool = Pool(meta.ncpu)
        queue_size = meta.ncpu
    else:
        meta.ncpu = 1
        pool = None
        queue_size = None

    print('Run dynesty...')
    if meta.run_dynamic:
        sampler = dynesty.DynamicNestedSampler(loglike, ptform, ndim, pool=pool, queue_size=queue_size,
                                               logl_args = l_args, ptform_args = p_args,
                                               update_interval=float(ndim), bound=meta.run_bound,
                                               sample=meta.run_sample)
        sampler.run_nested(dlogz_init=meta.run_dlogz_init, nlive_init=meta.run_nlive_init,
                           nlive_batch=meta.run_nlive_batch, maxbatch=meta.run_maxbatch)
    else:
        sampler = dynesty.NestedSampler(loglike, ptform, ndim, pool=pool, queue_size=queue_size,
                                        logl_args = l_args, ptform_args = p_args,
                                        update_interval=float(ndim), nlive=meta.run_nlive, bound=meta.run_bound,
                                        sample=meta.run_sample)
        sampler.run_nested(dlogz=meta.run_dlogz, print_progress=True)

    results = sampler.results

    # Closing multiprocessing
    if meta.ncpu > 1:
        pool.close()
        pool.join()

    # Dump the samples into a file using pickle
    with open(meta.workdir / meta.fitdir / 'nested_res' /
              f'nested_out_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.p', "wb") as pickle_file:
        pickle.dump(results, pickle_file)
    #results.summary()

    labels = meta.labels
    samples, weights = results.samples, np.exp(results.logwt - results.logz[-1])
    mean, cov = dyfunc.mean_and_cov(samples, weights)
    new_samples = dyfunc.resample_equal(samples, weights)

    # Saving plots
    plots.dyplot_runplot(results, meta)
    plots.dyplot_traceplot(results, meta)
    plots.dyplot_cornerplot(results, meta)

    # Determine median and 16th and 84th percentiles
    medians = []
    errors_lower = []
    errors_upper = []
    for i in range(ndim):
        q = util.quantile(new_samples[:, i], [0.16, 0.5, 0.84])
        medians.append(q[1])
        errors_lower.append(abs(q[1] - q[0]))
        errors_upper.append(abs(q[2] - q[1]))

    # Saving sampling results into txt files
    f_mcmc = open(meta.workdir / meta.fitdir / 'nested_res' /
                  f"nested_res_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.txt", 'w')
    for row in zip(errors_lower, medians, errors_upper, labels):
        print('{0: >8}: '.format(row[3]), '{0: >24} '.format(row[1]),
              '{0: >24} '.format(row[0]), '{0: >24} '.format(row[2]), file=f_mcmc)
    f_mcmc.close()

    updated_params = util.format_params_for_Model(medians, params, nvisit, fixed_array, tied_array, free_array)
    fit = model.fit(data, updated_params)
    util.append_fit_output(fit, meta, fitter='nested', medians=medians)

    # Saving plots
    plots.plot_fit_lc2(data, fit, meta, nested=True)
    plots.rmsplot(model, data, meta, fitter='nested')

    if meta.s30_fit_white:
        outfile = open(meta.workdir / meta.fitdir / 'white_systematics_nested.txt', "w")
        for i in range(len(fit.all_sys)):
            print(fit.all_sys[i], file=outfile)
        print('Saved white_systematics.txt file for nested sampling run')
        outfile.close()
    return medians, errors_lower, errors_upper, fit


def ptform(u, data):
    """Transforms the priors, which is needed for dynesty. """
    p = np.zeros_like(u)
    n = len(data.prior)
    for i in range(n):
        if data.prior[i][0] == 'U':
            p[i] = transform_uniform(u[i], data.prior[i][1], data.prior[i][2])
        if data.prior[i][0] == 'N':
            p[i] = transform_normal(u[i], data.prior[i][1], data.prior[i][2])
    return p


def loglike(x, params, data, model, nvisit,
            fixed_array, tied_array, free_array):
    """Calculates the log-likelihood."""
    updated_params = util.format_params_for_Model(x, params, nvisit, fixed_array, tied_array, free_array)
    if 'uncmulti' in data.s30_myfuncs:
        data.err = updated_params[-1] * data.err_notrescaled
    fit = model.fit(data, updated_params)
    return fit.ln_like
