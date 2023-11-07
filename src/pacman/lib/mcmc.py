import numpy as np
import pickle
import emcee
from multiprocessing import Pool
from . import plots
from . import util
from . import read_fit_par


def mcmc_fit(data, model, params, file_name, meta, fit_par):
    """
    Calls the emcee package and does the sampling.
    """
    nvisit = int(meta.nvisit)

    # Create the mcmc_res directory
    util.create_res_dir(meta)

    # Setting up parameters for sampler
    theta = util.format_params_for_sampling(params, meta, fit_par)
    ndim, nwalkers = len(theta), meta.run_nwalkers

    fixed_array = np.array(fit_par['fixed'])
    tied_array = np.array(fit_par['tied'])
    free_array = util.return_free_array(nvisit, fixed_array, tied_array)
    l_args = [params, data, model, nvisit, fixed_array, tied_array, free_array]

    # Setting up multiprocessing
    if hasattr(meta, 'ncpu') and meta.ncpu > 1:
        print('Using multiprocessing...')
        pool = Pool(meta.ncpu)
    else:
        meta.ncpu = 1
        pool = None

    print('Run emcee...')
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = l_args, pool=pool)
    step_size = read_fit_par.get_step_size(data, params, meta, fit_par)
    pos = [theta + np.array(step_size)*np.random.randn(ndim) for i in range(nwalkers)]
    pos = np.array(pos)

    sampler.run_mcmc(pos, meta.run_nsteps, progress=True)

    # Dump the samples into a file using pickle
    with open(meta.workdir + meta.fitdir + '/mcmc_res/' +
              '/mcmc_out_bin{0}_wvl{1:0.3f}.p'.format(meta.s30_file_counter, meta.wavelength), "wb") as pickle_file:
        pickle.dump([data, params, sampler.chain], pickle_file)
    nburn = meta.run_nburn

    labels = meta.labels

    # Thin out the Samples if a huge amount was used.
    # This was introduced as the run would break when saving the plot due to its big size
    if meta.run_nsteps * meta.run_nwalkers > 1000000:
        thin_corner = int((meta.run_nsteps - meta.run_nburn) * meta.run_nwalkers // 100000)
        print('Note: Big Corner plot with many steps. Thinning Plot by factor: {0}'.format(thin_corner))
    else:
        thin_corner = 1

    samples = sampler.chain[:, nburn::thin_corner, :].reshape((-1, ndim))

    # Closing multiprocessing
    if meta.ncpu > 1:
        pool.close()
        pool.join()

    # Saving plots
    plots.mcmc_pairs(samples, params, meta, fit_par, data)
    plots.mcmc_chains(ndim, sampler, 0, labels, meta)
    plots.mcmc_chains(ndim, sampler, nburn, labels, meta)

    # Determine median and 16th and 84th percentiles
    medians = []
    errors_lower = []
    errors_upper = []
    for i in range(len(theta)):
        q = util.quantile(samples[:, i], [0.16, 0.5, 0.84])
        medians.append(q[1])
        errors_lower.append(abs(q[1] - q[0]))
        errors_upper.append(abs(q[2] - q[1]))

    # Saving sampling results into txt files
    f_mcmc = open(meta.workdir + meta.fitdir + '/mcmc_res/' +
                  "/mcmc_res_bin{0}_wvl{1:0.3f}.txt".format(meta.s30_file_counter, meta.wavelength), 'w')
    for row in zip(errors_lower, medians, errors_upper, labels):
        print(f'{row[3]: >8}: {row[1]: >24} {row[0]: >24} {row[2]: >24} ', file=f_mcmc)
    f_mcmc.close()

    updated_params = util.format_params_for_Model(medians, params, nvisit, fixed_array, tied_array, free_array)
    fit = model.fit(data, updated_params)
    util.append_fit_output(fit, meta, fitter='mcmc', medians=medians)

    # Saving plots
    plots.plot_fit_lc2(data, fit, meta, mcmc=True)
    plots.rmsplot(model, data, meta, fitter='mcmc')

    if meta.s30_fit_white:
        outfile = open(meta.workdir + meta.fitdir + '/white_systematics_mcmc.txt', "w")
        for i in range(len(fit.all_sys)): print(fit.all_sys[i], file=outfile)
        print('Saved white_systematics.txt file for mcmc run')
        outfile.close()

    #samples_auto1 = sampler.get_chain(discard=nburn, flat=True)
    #tau = emcee.autocorr.integrated_time(samples_auto1, quiet=True)
    #print(f"Autocorrelation time: {tau}")

    #tau = sampler.get_autocorr_time(discard=nburn)
    #print(f"Autocorrelation time: {tau}")

    # Print the Autocorrelation Time
    tau = sampler.get_autocorr_time(quiet=True)
    print("Autocorrelation time: ", tau)

    return medians, errors_lower, errors_upper, fit


def lnprior(theta, data):
    """
    Calculate the log-prior.
    """
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


def lnprob(theta, params, data, model, nvisit, fixed_array, tied_array, free_array):
    """
    Calculates the log-likelihood.
    """
    updated_params = util.format_params_for_Model(theta, params, nvisit, fixed_array, tied_array, free_array)
    if 'uncmulti' in data.s30_myfuncs:
        data.err = updated_params[-1] * data.err_notrescaled
    lp = lnprior(theta, data)
    if lp == -np.inf:  # if the likelihood from the priors is already -inf, dont evaluate the function
        return lp
    fit = model.fit(data, updated_params)
    return fit.ln_like + lp


#ORDER
#mcmc_fit
#format_params_for_mcmc
#mcmc_fit
#lnprob
#format_params_for_Model
#lnprob
#lnprior
#lnprob
