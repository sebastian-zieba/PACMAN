import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')

import wfc3_reduction.reduction.s02_barycorr as s02
from wfc3_reduction.lib.update_meta import update_meta

eventlabel = 'KELT11_Hubble15926'
workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-12-30_17-49-34_KELT11_Hubble15926/'

update_meta(eventlabel, workdir)
meta = s02.run02(eventlabel, workdir)

