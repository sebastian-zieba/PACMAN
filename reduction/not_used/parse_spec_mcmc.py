import sys                                                                      
sys.path.insert(0, './fit_funcs/')                                              
sys.path.insert(0, 'fit_funcs/models')                                          
import pickle                                                                   
import os                                                                       
import glob                                                                     
import numpy as np                                                              
from read_data import Data                                                      
from model import Model                                                         
import matplotlib.pyplot as plt 

def quantile(x, q): return np.percentile(x, [100. * qi for qi in q]) 

path = "./spec_mcmc/"

mcmc = glob.glob(os.path.join(path, "mcmc*.p"))
lsq = glob.glob(os.path.join(path, "lsq*.p"))

print "FIXME setting ndim by hand"
print "AHHHHHHHHHHH"

ndim = 11

waves, ys, yes = [], [], []

for m, l in zip(mcmc, lsq):
    data, params, chain = pickle.load(open(m, "r"))
    data, model = pickle.load(open(l, "r"))


    samples = chain[:, 1000:, :].reshape((-1, ndim))
    
    medians, errors = [], []

    for i in range(ndim):
            q = quantile(samples[:, i], [0.16, 0.5, 0.84])
            medians.append(q[1])
            errors.append(q[2] - q[1])
    
    y, ye = medians[0], errors[0]
    waves.append(data.wavelength)
    ys.append(y)
    yes.append(ye)
    print data.wavelength,  y, ye 

plt.errorbar(np.array(waves), np.array(ys), yerr = np.array(yes), fmt = '.k')
plt.ylim(0, 2.e-3)
plt.show()



