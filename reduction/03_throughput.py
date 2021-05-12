import numpy as np
import matplotlib.pyplot as plt
import pysynphot as S 

bp_g102 = S.ObsBandpass('wfc3,ir,g102') #FIXME add mjd check to download time specific bandpass
bp_g141 = S.ObsBandpass('wfc3,ir,g141') 

wvl_g102, thr_g102 = bp_g102.binset, bp_g102(bp_g102.binset)
wvl_g141, thr_g141 = bp_g141.binset, bp_g141(bp_g141.binset)

np.savetxt('./ancil/throughputs_and_spectra/g102_throughput.txt', list(zip(wvl_g102, thr_g102)))
np.savetxt('./ancil/throughputs_and_spectra/g141_throughput.txt', list(zip(wvl_g141, thr_g141)))

plt.plot(wvl_g102, thr_g102) 
plt.plot(wvl_g141, thr_g141) 

plt.show()

