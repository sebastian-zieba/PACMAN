import sys
sys.path.insert(0,'..')
import numpy as np

def exponential_visit_full(t, data, params, visit = 0):
    scale, C, exp1, exp2 = params
    scale = scale[visit]
    C = C[visit]
    C = 10.**C
    exp1 = exp1[visit]
    exp2 = exp2[visit]

    t_vis = data.t_vis[data.vis_idx[visit]]

    return (1. + scale*data.scan_direction[data.vis_idx[visit]]) * C + np.exp(- t_vis * exp1) * exp2
