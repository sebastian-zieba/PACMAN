import sys
sys.path.insert(0,'..')
import celerite
from celerite import terms
import numpy as np
import matplotlib.pyplot as plt


def gp_matern32(idx, data, resid, params):
    log_rho, log_sigma = params
    err = data.err/data.flux

    t = data.t_vis[idx]/60/60
    kernel = terms.Matern32Term(log_rho=log_rho, log_sigma=log_sigma)

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
