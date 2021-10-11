import sys, os, time
#sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
from importlib import reload

import wfc3_reduction.reduction.s00_table as s00
import wfc3_reduction.reduction.s01_horizons as s01
import wfc3_reduction.reduction.s02_barycorr as s02
import wfc3_reduction.reduction.s03_refspectra as s03
import wfc3_reduction.reduction.s10_direct_images as s10
import wfc3_reduction.reduction.s20_extract as s20
import wfc3_reduction.reduction.s21_ld_inputmaker as s21
import wfc3_reduction.reduction.s22_bin_spectroscopic_lc as s22

eventlabel = 'L-98-59_Hubble15856'
#eventlabel = 'KELT11_Hubble15926'

reload(s00)
meta = s00.run00(eventlabel)

reload(s01)
meta = s01.run01(meta.eventlabel, meta.workdir, meta=meta)

reload(s02)
meta = s02.run02(meta.eventlabel, meta.workdir, meta=meta)

reload(s03)
meta = s03.run03(meta.eventlabel, meta.workdir, meta=meta)

reload(s10)
meta = s10.run10(meta.eventlabel, meta.workdir, meta=meta)

reload(s20)
meta = s20.run20(meta.eventlabel, meta.workdir, meta=meta)

#reload(s21)
#meta = s21.run21(meta.eventlabel, meta.workdir, meta=meta)

#reload(s22)
#meta = s22.run22(meta.eventlabel, meta.workdir, meta=meta)
