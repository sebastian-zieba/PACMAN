import sys
sys.path.insert(0,'..')
import numpy as np


def upstream_downstream(t, data, params, visit = 0):
    scale = params
    scale = scale[0][visit]

    return 1. + scale*data.scan_direction[data.vis_idx[visit]]
