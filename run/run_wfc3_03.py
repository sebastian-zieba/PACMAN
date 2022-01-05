import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')

import wfc3_reduction.reduction.s03_refspectra as s03
from wfc3_reduction.lib.update_meta import update_meta

eventlabel = 'KELT11_Hubble15926'
workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-12-30_17-49-34_KELT11_Hubble15926/'

update_meta(eventlabel, workdir)
meta = s03.run03(eventlabel, workdir)

