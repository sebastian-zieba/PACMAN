import sys
import numpy as np

sys.path.insert(0,'..')


def sine1(t, data, params):
    a1, omega1, phi1  = params
    # FIXME: data.t_vis won't work if there are multiple visits
    return 1. + a1*np.sin(omega1*data.t_vis + phi1)
