import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload

import wfc3_reduction.reduction.s20_extract as s20
from wfc3_reduction.lib.reload_meta import reload_meta
eventlabel = 'L-98-59_Hubble15856'
workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-13_16-05-26_L-98-59_Hubble15856/'
reload_meta(workdir, eventlabel)
reload(s20)
s20.run20(eventlabel, workdir)
