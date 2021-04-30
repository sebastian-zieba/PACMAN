import sys
sys.path.insert(0,'..')
from read_data import Data
import numpy as np

def sine1(t, data, params):
    a1, omega1, phi1  = params
    #FIXME data.t_vis won't work if there are multiple visits
    return 1. + a1*np.sin(omega1*data.t_vis + phi1)
