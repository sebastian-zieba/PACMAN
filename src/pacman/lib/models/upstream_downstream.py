import sys

sys.path.insert(0, '..')


def upstream_downstream(t, data, params, visit: float = 0.):
    scale = params
    scale = scale[0][visit]
    return 1. + scale*data.scan_direction[data.vis_idx[visit]]
