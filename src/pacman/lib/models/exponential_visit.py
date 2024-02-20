import sys
import numpy as np

sys.path.insert(0, '..')

def exponential_visit(t, data, params, visit: float = 0.):
    exp1, exp2 = params
    exp1 = exp1[visit]
    exp2 = exp2[visit]

    idx = data.vis_idx[visit]
    t_vis = data.t_vis[data.vis_idx[visit]]

    return (1. - exp1 * np.exp(-t_vis/exp2))
    #return (1. - np.exp(-exp1 * t_vis - exp2)) # alternative parametrization
