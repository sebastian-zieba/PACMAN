import sys
sys.path.insert(0,'..')
import numpy as np
from read_data import Data

def constant(t, data, params, visit = 0):
    C = params
    C = C[0][visit]
    C = 10.**C

    return 1. + C*np.ones_like(t)
