import sys

sys.path.insert(0,'..')


def polynomial1(t, data, params, visit: float = 0.):
    v = params
    v = v[0][visit]

    t_vis = data.t_vis[data.vis_idx[visit]]

    return 1. + v*t_vis
