import numpy as np
import pickle
from datetime import datetime
from scipy.stats import norm
import dynesty
import inspect
from dynesty import plotting as dyplot
import matplotlib.pyplot as plt

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

def format_params_for_mcmc(params, obs_par, fit_par):	#FIXME: make sure this works for cases when nvisit>1
    nvisit = int(obs_par['nvisit'])				
    theta = []

    for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "false":
                    if fit_par['tied'][i].lower() == "true": theta.append(params[i*nvisit])
                    else: 
                            for j in range(nvisit): theta.append(params[i*nvisit+j])
    return np.array(theta)


def mcmc_output(samples, params, obs_par, fit_par, data):	#FIXME: make sure this works for cases when nvisit>1
    nvisit = int(obs_par['nvisit'])				
    labels = []

    for i in range(len(fit_par)):
            if fit_par['fixed'][i].lower() == "false":
                    if fit_par['tied'][i].lower() == "true": labels.append(fit_par['parameter'][i])
                    else: 
                            for j in range(nvisit): labels.append(fit_par['parameter'][i]+str(j))
    fig = corner.corner(samples, labels=labels, show_titles=True)
    current_time = datetime.now().time()
    figname = "pairs_"+current_time.isoformat()+".png"
    fig.savefig(figname)


def format_params_for_Model(theta, params, obs_par, fit_par):
    nvisit = int(obs_par['nvisit'])
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

def nested_sample(data, model, params, file_name, obs_par, fit_par):
    x = format_params_for_mcmc(params, obs_par, fit_par)	

    ndim = len(x) 

    l_args = [params, data, model, obs_par, fit_par]
    p_args = [data]
    
    dsampler = dynesty.DynamicNestedSampler(loglike, ptform, ndim, 
                                            logl_args = l_args,
                                            ptform_args = p_args,
                                            update_interval=float(ndim))
    dsampler.run_nested(wt_kwargs={'pfrac': 1.0})#, maxiter = 2000)

    results = dsampler.results

    pickle.dump(results, open( "nested_results.p", "wb" ) )
    #results.summary()

    # Plot a summary of the run.
    """rfig, raxes = dyplot.runplot(results)

    # Plot traces and 1-D marginalized posteriors.
    tfig, taxes = dyplot.traceplot(results)

    # Plot the 2-D marginalized posteriors.
    cfig, caxes = dyplot.cornerplot(results)

    plt.show()"""


    return 0.


def ptform(u, data):
    p = np.zeros_like(u) 
    n = len(data.prior)
    for i in range(n):
        if data.prior[i][0] == 'U':  p[i] = transform_uniform(u[i], 
                                            data.prior[i][1],data.prior[i][2])
        if data.prior[i][0] == 'N':  p[i] = transform_normal(u[i], 
                                            data.prior[i][1],data.prior[i][2])
    return p
    


def loglike(x, params, data, model, obs_par, fit_par):
    updated_params = format_params_for_Model(x, params, obs_par, fit_par)
    fit = model.fit(data, updated_params)
    return fit.ln_like 
