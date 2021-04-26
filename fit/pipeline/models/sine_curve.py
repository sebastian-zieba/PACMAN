import numpy as np

def get_phaselc(t, p, data, v_num):
	return 1.+p.amp1[v_num]*np.cos(2.*np.pi*(t-p.theta1[v_num])/p.per[v_num]) + p.amp2[v_num]*np.cos(4.*np.pi*(t-p.theta2[v_num])/p.per[v_num])

