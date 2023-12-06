import sys
import numpy as np

sys.path.insert(0,'..')


def model_ramp(t, data, params, visit: float = 0.):
    r1, r2, r3 = params
    r1 = r1[visit]
    r2 = r2[visit]
    r3 = r3[visit]

    t_orb = data.t_orb[data.vis_idx[visit]]
    t_delay = data.t_delay[data.vis_idx[visit]]

    return 1.0 - np.exp(-r1*t_orb - r2 - r3*t_delay)
