import sys
import numpy as np


sys.path.insert(0,'..')


def constant(t, data, params, visit = 0):
    C = params
    C = C[0][visit]
    C = 10.**C
    return 1e-10 + C*np.ones_like(t)
