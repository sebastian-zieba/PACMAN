import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload

import wfc3_reduction.reduction.s21_ld_inputmaker as s21

eventlabel = 'L-98-59_Hubble15856'

workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-03_12-22-13_L-98-59_Hubble15856/' #[0]
#workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-04_03-19-10_L-98-59_Hubble15856/' #[0,2]
#workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-08-04_04-00-53_L-98-59_Hubble15856/' #[all]

s21.run21(eventlabel, workdir)