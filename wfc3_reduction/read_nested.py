import dynesty
import pickle
import numpy as np
import corner
import matplotlib.pyplot as plt
from dynesty import utils as dyfunc
from pylab import  hist
from dynesty import plotting as dyplot
from dynesty.utils import quantile


d= pickle.load(open("nested_results.p", "rb"))

res2 = d

samples, weights = res2.samples, np.exp(res2.logwt - res2.logz[-1])
mean, cov = dyfunc.mean_and_cov(samples, weights)

new_samples = dyfunc.resample_equal(samples, weights)

#print new_samples[:,0], new_samples[:,0].min(), new_samples[:,0].max()
#hist(new_samples[:,0])
#plt.show()

#corner.corner(new_samples[::10])

new_samples = dyfunc.resample_equal(samples, weights)
q = quantile(new_samples[:, 0], [0.16, 0.5, 0.84])
#2458453.4 is toffset from obs_par.txt
print("t0 = ", 2458453.4 + q[1], "+", q[2] - q[1], "-", q[1] - q[0])

q = quantile(new_samples[:, 1], [0.16, 0.5, 0.84])
print("a/Rs = ", q[1], q[1] - q[0])

q = quantile(new_samples[:, 2], [0.16, 0.5, 0.84])
print("inc = ", q[1], q[1] - q[0])

corner.corner(new_samples)

#corner.corner(samples)

plt.savefig("corners.png")
