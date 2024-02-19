import sys
sys.path.insert(0,'..')
import numpy as np

def logarithmic_visit_full(t, data, params, visit = 0):
    scale, C, log1, log2 = params
    scale = scale[visit]
    C = C[visit]
    C = 10.**C
    log1 = log1[visit]
    log2 = log2[visit]

    t_vis = data.t_vis[data.vis_idx[visit]]

    return (1. + scale*data.scan_direction[data.vis_idx[visit]]) * C + np.log(t_vis + 1/log1) * log2
