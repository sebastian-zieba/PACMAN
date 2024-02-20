import sys
sys.path.insert(0,'..')
import numpy as np

def logarithmic_visit(t, data, params, visit = 0):
    log1, log2 = params
    log1 = log1[visit]
    log2 = log2[visit]

    idx = data.vis_idx[visit]
    t_vis = data.t_vis[data.vis_idx[visit]]

    return (1. - log1 * np.log(t_vis + log2))
    #return (1. - log1 * np.log(t_vis + 10**log2)) # alternative parametrization
