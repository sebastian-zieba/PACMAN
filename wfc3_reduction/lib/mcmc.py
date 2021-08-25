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
    #print('theta', theta)
    return np.array(theta)


def get_step_size(params, meta, fit_par):
    nvisit = int(meta.nvisit)
    step_size = []

    if meta.fit_par_new == False:
        for i in range(len(fit_par)):
                if fit_par['fixed'][i].lower() == "false":
                        if fit_par['tied'][i].lower() == "true": step_size.append(fit_par['step_size'][i])
                        else:
                                for j in range(nvisit): step_size.append(fit_par['step_size'][i])
    else:
        ii = 0
        for i in range(int(len(params)/nvisit)):
                if fit_par['fixed'][ii].lower() == "false":
                        if str(fit_par['tied'][ii]) == "-1":
                            step_size.append(fit_par['step_size'][ii])
                            ii = ii + 1
                        else:
                            for j in range(nvisit):
                                step_size.append(fit_par['step_size'][ii])
                                ii = ii + 1
                else:
                    ii = ii + 1

    return np.array(step_size)

def labels_gen(params, meta, fit_par):
    nvisit = int(meta.nvisit)
    labels = []

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
    labels = labels_gen(params, meta, fit_par)

    fig = corner.corner(samples, labels=labels, show_titles=True,quantiles=[0.16, 0.5, 0.84],title_fmt='.4')

    figname = meta.workdir + meta.fitdir + "/pairs_{0}_".format(meta.run_file.split('/')[-1]) + meta.fittime + ".png"
    fig.savefig(figname)


def format_params_for_Model(theta, params, meta, fit_par):
    nvisit = int(meta.nvisit)
    params_updated = []
    iter = 0									#this should be a more informative name FIXME

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

   # print('params', params)
    #print('params_updated', params_updated)
    #print(len(params))
    #print(len(params_updated))
    return np.array(params_updated)

def mcmc_fit(data, model, params, file_name, meta, fit_par):
    theta = format_params_for_mcmc(params, meta, fit_par)

    ndim, nwalkers = len(theta), meta.run_nwalkers
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (params, data, model, meta, fit_par))


    step_size = get_step_size(params, meta, fit_par)
    print(step_size)
    print(theta)
    [print("({0},{1})".format(i, j)) for i, j in zip(theta,step_size)]
    pos = [theta + np.array(step_size)*np.random.randn(ndim) for i in range(nwalkers)]


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
    labels = labels_gen(params, meta, fit_par)

    plots.plot_chains(ndim, sampler, 0, labels, meta)
    plots.plot_chains(ndim, sampler, nburn, labels, meta)


    medians = []
    errors_lower = []
    errors_upper = []
    errors_mean = []
    for i in range(len(theta)):
            q = quantile(samples[:, i], [0.16, 0.5, 0.84])
            medians.append(q[1])
            errors_lower.append(abs(q[1] - q[0]))
            errors_upper.append(abs(q[2] - q[1]))
            errors_mean.append(abs((q[2] - q[0])/2))

    f_mcmc = open(meta.workdir + meta.fitdir + "/mcmc_res_{0}_{1}.txt".format(meta.run_file.split('/')[-1], meta.fittime), 'w')
    for row in zip(errors_lower, medians, errors_upper, labels):
        #print("{: <8} +/-  {: <8}".format(*row))
        xl = ufloat(row[1], row[0])
        xu = ufloat(row[1], row[2])
        print('{0: >12}: '.format(row[3]), '{0: >12} '.format(row[0]), '{0: >12} '.format(row[1]), '{0: >12} '.format(row[2]), file=f_mcmc)
        print('2:', file=f_mcmc)
        print('lower', file=f_mcmc)
        print('{:13.2uS}'.format(xl), file=f_mcmc)
        print('{:13.2u}'.format(xl), file=f_mcmc)
        print('upper', file=f_mcmc)
        print('{:13.2uS}'.format(xu), file=f_mcmc)
        print('{:13.2u}'.format(xu), file=f_mcmc)

        print('3:', file=f_mcmc)
        print('lower', file=f_mcmc)
        print('{:13.3uS}'.format(xl), file=f_mcmc)
        print('{:13.3u}'.format(xl), file=f_mcmc)
        print('upper', file=f_mcmc)
        print('{:13.3uS}'.format(xu), file=f_mcmc)
        print('{:13.3u}'.format(xu), file=f_mcmc)
    f_mcmc.close()

    updated_params = format_params_for_Model(medians, params, meta, fit_par)
    fit = model.fit(data, updated_params)
    print(fit.rms)

    plots.plot_fit_lc2(data, fit, meta, mcmc=True)

    return data.wavelength, medians[0], errors_mean,samples


def lnprior(theta, data):
    lnprior_prob = 0.
    n = len(data.prior)
    for i in range(n):
        #print(lnprior_prob)
        if data.prior[i][0] == 'U': 
            if np.logical_or(theta[i] < data.prior[i][1], 
              theta[i] > data.prior[i][2]): lnprior_prob += - np.inf
        if data.prior[i][0] == 'N': 
            lnprior_prob -= 0.5*(np.sum(((theta[i] - 
              data.prior[i][1])/data.prior[i][2])**2 + 
              np.log(2.0*np.pi*(data.prior[i][2])**2)))
    #print('lnprior_prob', lnprior_prob)
    return lnprior_prob
    


def lnprob(theta, params, data, model, meta, fit_par):
    updated_params = format_params_for_Model(theta, params, meta, fit_par)
    #print('updated_params', updated_params[12*5:14*5])
    fit = model.fit(data, updated_params)
    #print('fit', fit)
    lp = lnprior(theta, data)
    #print('lp', lp)
    return fit.ln_like + lp



#mcmc_fit
#format_params_for_mcmc
#mcmc_fit
#lnprob
#format_params_for_Model
#lnprob
#lnprior
#lnprob



