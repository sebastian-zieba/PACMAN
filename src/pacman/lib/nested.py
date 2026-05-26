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
    untied_array = util.return_untied_array(nvisit, fixed_array, tied_array)
    l_args = [params, data, model, nvisit, fixed_array, tied_array, free_array, untied_array]
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
    with (meta.workdir / meta.fitdir / 'nested_res' /
          f'nested_out_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.p')\
            .open("wb") as pickle_file:
        pickle.dump(results, pickle_file)
    #results.summary()

    labels = meta.labels
    samples, weights = results.samples, np.exp(results.logwt - results.logz[-1])
    mean, cov = dyfunc.mean_and_cov(samples, weights)
    new_samples = dyfunc.resample_equal(samples, weights)

    best_idx_ml = np.argmax(results.logl)
    best_sample_ml = results.samples[best_idx_ml]

    # Determine median and 16th and 84th percentiles
    p16_list = []
    p50_list = []
    p84_list = []
    errors_lower = []
    errors_upper = []
    for i in range(ndim):
        q = util.quantile(new_samples[:, i], [0.16, 0.5, 0.84])
        p16, p50, p84 = q
        p16_list.append(p16)
        p50_list.append(p50)
        p84_list.append(p84)
        errors_lower.append(p50 - p16)
        errors_upper.append(p84 - p50)
    medians = p50_list

    plots.nested_pairs(
        new_samples,
        meta.labels,
        meta,
        median_vals=medians,
        ml_vals=best_sample_ml,
        n_sampled=ndim,
    )

    # Saving plots
    plots.dyplot_runplot(results, meta)
    plots.dyplot_traceplot(results, meta)
    plots.dyplot_cornerplot(results, meta)


    # Saving sampling results into txt files
    with (meta.workdir / meta.fitdir / 'nested_res' /
        f"nested_res_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.txt").open('w', encoding='utf-8') as f_nested:
        print(
            f"{'parameter':<20} {'p50':>14} {'p16':>14} {'p84':>14} {'minus':>14} {'plus':>14} {'maxL':>14}",
            file=f_nested
        )
        for label, p50, p16, p84, minus, plus, ml in zip(
            labels, p50_list, p16_list, p84_list, errors_lower, errors_upper, best_sample_ml
        ):
            print(
                f"{label:<20} {p50:14.7g} {p16:14.7g} {p84:14.7g} {minus:14.7g} {plus:14.7g} {ml:14.7g}",
                file=f_nested
            )

    updated_params = util.format_params_for_Model(medians, params, nvisit, fixed_array, tied_array, free_array, untied_array)
    if "uncmulti" in data.s30_myfuncs:
        util.apply_uncmulti(data, updated_params)
    fit = model.fit(data, updated_params)
    util.append_fit_output(fit, meta, fitter='nested', medians=medians)

    # Saving plots
    plots.plot_fit_lc2(data, fit, meta, nested=True)
    plots.rmsplot(model, data, meta, fitter='nested')
    plots.save_astrolc_data(data, model, meta, fitter='nested')

    if meta.s30_fit_white:
        with (meta.workdir / meta.fitdir / 'white_systematics_nested.txt')\
                .open("w", encoding='utf-8') as outfile:
            for i in range(len(fit.all_sys)):
                print(fit.all_sys[i], file=outfile)
        print('Saved white_systematics.txt file for nested sampling run')
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
            fixed_array, tied_array, free_array, untied_array):
    """Calculates the log-likelihood."""
    updated_params = util.format_params_for_Model(x, params, nvisit, fixed_array, tied_array, free_array, untied_array)
    if "uncmulti" in data.s30_myfuncs:
        util.apply_uncmulti(data, updated_params)
    fit = model.fit(data, updated_params)
    return fit.ln_like
