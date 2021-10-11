import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload

import wfc3_reduction.reduction.s20_extract as s20
from wfc3_reduction.lib.reload_meta import reload_meta

#eventlabel = 'L-98-59_Hubble15856'
eventlabel = 'KELT11_Hubble15926'

#workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-21_02-00-29_L-98-59_Hubble15856_spec'
#workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-09-08_09-38-00_L-98-59_Hubble15856/'

workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-09-05_22-39-58_KELT11_Hubble15926_spec2/'

reload_meta(eventlabel, workdir)
reload(s20)
s20.run20(eventlabel, workdir)
