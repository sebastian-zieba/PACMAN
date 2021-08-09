import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload

import wfc3_reduction.reduction.s20_extract as s20

eventlabel = 'L-98-59_Hubble15856'

workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-07-28_23-40-22_L-98-59_Hubble15856'

reload(s20)
s20.run20(eventlabel, workdir)
