import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')

import wfc3_reduction.reduction.s00_table as s00

eventlabel = 'KELT11_Hubble15926'

meta = s00.run00(eventlabel)

