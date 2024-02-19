import sys
sys.path.insert(0,'..')
import numpy as np

def polynomial1_full(t, data, params, visit = 0):
    scale, C, v = params
    scale = scale[visit]
    C = C[visit]
    C = 10.**C
    v = v[visit]

    t_vis = data.t_vis[data.vis_idx[visit]]

    return (1. + scale*data.scan_direction[data.vis_idx[visit]]) * C + v*t_vis
