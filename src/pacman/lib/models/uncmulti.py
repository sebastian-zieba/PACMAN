import numpy as np


def uncmulti(t, data, params, visit = 0):
    val = params
    val = val[0][visit]
    #data.err = data.err * val
    return np.ones_like(t)
