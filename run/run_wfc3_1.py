import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload

import wfc3_reduction.reduction.s10_direct_images as s10

eventlabel = 'L-98-59_Hubble15856'

workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-07-27_23-43-18_L-98-59_Hubble15856'

reload(s10)
meta = s10.run10(eventlabel, workdir)

