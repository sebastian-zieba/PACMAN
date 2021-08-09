import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload

import wfc3_reduction.reduction.s00_table as s00

eventlabel = 'L-98-59_Hubble15856'

reload(s00)
meta = s00.run00(eventlabel)

