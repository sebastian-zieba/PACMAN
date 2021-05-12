import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.optimize import leastsq


def residuals(params, x1, y1, x2, y2):
    a, b = params    
    x1=np.array(x1)
    x2=np.array(x2)
    y1=np.array(y1)
    y2=np.array(y2)

    f = interp1d(x1, y1, kind='cubic')

    fit = f(a+b*x2)

    return fit - y2

def rms(x):
    return sum(x**2)


cmin, cmax = 270, 451

data=np.loadtxt('/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/reduction/not_used/L89-59-exp1.txt').T
model=np.loadtxt('/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/reduction/ancil/throughputs_and_spectra/thoughput_times_spectrum.txt').T


datax = data[0][cmin:cmax]
datay = data[1][cmin:cmax]/max(data[1][cmin:cmax])

modelx=np.concatenate((np.linspace(0.001, min(model[0]), 100, endpoint=False), model[0], np.linspace(max(model[0])+0.01, 100, 100, endpoint=False)))
modely=np.concatenate((np.zeros(100), model[1]/max(model[1]),np.zeros(100)))



p0=[-0.4, 1/200]

leastsq_res = leastsq(residuals, p0, args=(modelx, modely, datax, datay))[0]

print(leastsq_res)


#a_list=np.linspace(-0.6,0, 30)
#b_list=np.linspace(1/250, 1/150, 30)
#res = [] 
#combos = np.array(np.meshgrid(a_list, b_list)).T.reshape(-1, 2)

#for combo in combos:    
#    res.append(rms(residuals(combo, model0, model1, data[0][280:440], data[1][280:440]/max(data[1]))))

#print(combos[np.argmin(res)][0])
#print(combos[np.argmin(res)][1])


plt.plot(datax*leastsq_res[1]+leastsq_res[0], datay, label = 'spectrum fit, wvl = {0:.5g}+{1:.5g}*pixel'.format(leastsq_res[0],leastsq_res[1]))
plt.plot(modelx, modely, label = 'throughput * spectrum')
plt.xlim(0,2)
plt.legend()
plt.savefig('comp1.png')
#plt.show()


