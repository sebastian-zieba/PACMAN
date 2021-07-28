import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload

import wfc3_reduction.reduction.s30_run as s30

eventlabel = 'L-98-59_Hubble15856'

#workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-07-28_01-15-19_L-98-59_Hubble15856/'
workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-07-28_02-40-16_L-98-59_Hubble15856'

reload(s30)
s30.run30(eventlabel, workdir)
