import sys

import celerite
import numpy as np
from celerite import terms

sys.path.insert(0,'..')


def gp_sho(idx, data, resid, params):
    logQ, logw, logS, log_jit = params
    err = data.err/data.flux
    #logerr = np.log(np.median(err))
    logerr = log_jit

    #print "logQ, logw, logS", logQ, logw, logS

    mean_resid = np.mean(resid)
    #resid -= mean_resid

    t = data.t_vis[idx]/60/60

    kernel = (terms.SHOTerm(log_S0 = logS, log_Q = logQ, log_omega0 = logw) + 
              terms.JitterTerm(log_sigma = logerr)
              )

    gp = celerite.GP(kernel, fit_mean = True) 
    gp.compute(t, err, check_sorted = False)   #t must be ascending!


    mu = gp.predict(resid, t, return_cov = False)
    gp_lnlike = gp.log_likelihood(resid)

    # plt.errorbar(t, resid, err, fmt = '.k')
    # plt.plot(t, mu)
    #
    # x = np.linspace(np.min(t), np.max(t), 1000)
    # pred_mean, pred_var = gp.predict(resid, x, return_var = True)
    # pred_std = np.sqrt(pred_var)
    # plt.fill_between(x, pred_mean+pred_std, pred_mean-pred_std, color='blue', alpha=0.3)
    #
    # plt.show()

    #np.save("resid", [t,  resid, err])

    #return [1.0 + np.array(mu) + mean_resid, gp_lnlike] 
    return [1.0 + np.array(mu), gp_lnlike]
