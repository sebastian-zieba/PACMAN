import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def polynomial2(t, data, params, visit = 0):
    v, v2 = params
    v = v[visit]
    v2 = v2[visit]

    idx = data.vis_idx[visit]

    return 1. + v*data.t_vis[idx] + v2*(data.t_vis[idx])**2
