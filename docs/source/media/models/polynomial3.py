import sys
sys.path.insert(0,'..')
import numpy as np

def polynomial3(t, data, params, visit = 0):
    v, v2, v3 = params
    v = v[visit]
    v2 = v2[visit]
    v3 = v3[visit]

    idx = data.vis_idx[visit]

    return 1. + v*data.t_vis[idx] + v2*(data.t_vis[idx])**2 + v3*(data.t_vis[idx])**3
