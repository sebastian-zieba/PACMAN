import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload

import wfc3_reduction.reduction.s22_bin_spectroscopic_lc as s22

eventlabel = 'L-98-59_Hubble15856'

workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-07-22_06-58-22_L-98-59_Hubble15856'

s22.run22(eventlabel, workdir)