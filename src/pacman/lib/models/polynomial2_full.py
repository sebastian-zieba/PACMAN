import sys
sys.path.insert(0,'..')
import numpy as np

def polynomial2_full(t, data, params, visit = 0):
    scale, C, v, v2 = params
    scale = scale[visit]
    C = C[visit]
    C = 10.**C
    v = v[visit]
    v2 = v2[visit]

    t_vis = data.t_vis[data.vis_idx[visit]]

    return (1. + scale*data.scan_direction[data.vis_idx[visit]]) * C + v*t_vis + v2*(t_vis)**2
