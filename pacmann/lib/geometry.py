import numpy as np

#calculates the slope and intercept for the trace, given the position  of the direct image in physical pixels
#these coefficients are for the WFC3 G141 grism
def trace(X_ref, Y_ref):		
	BEAMA_i = 15 
	BEAMA_f = 196

	DYDX_0_0 = 1.96882E+00 
	DYDX_0_1 = 9.09159E-05
	DYDX_0_2 = -1.93260E-03

	DYDX_1_0 = 1.04275E-02 
	DYDX_1_1 = -7.96978E-06 
	DYDX_1_2 = -2.49607E-06  
	DYDX_1_3 = 1.45963E-09  
	DYDX_1_4 = 1.39757E-08  
	DYDX_1_5 = 4.84940E-10


	DYDX_0 = DYDX_0_0 + DYDX_0_1*X_ref + DYDX_0_2*Y_ref
	DYDX_1 = DYDX_1_0 + DYDX_1_1*X_ref + DYDX_1_2*Y_ref + DYDX_1_3*X_ref**2 + DYDX_1_4*X_ref*Y_ref + DYDX_1_5*Y_ref**2
	
	return [DYDX_0, DYDX_1]

#calculates coefficients for the dispersion solution
def dispersion(X_ref, Y_ref):		#X_ref and Y_ref are the centroid position in physical pixels
	DLDP_0_0 = 8.95431E+03   
	DLDP_0_1 = 9.35925E-02   
	DLDP_0_2 = 0.0

	DLDP_1_0 = 4.51423E+01
	DLDP_1_1 = 3.17239E-04  
	DLDP_1_2 = 2.17055E-03 
	DLDP_1_3 = -7.42504E-07 
	DLDP_1_4 = 3.48639E-07 
	DLDP_1_5 = 3.09213E-07

	DLDP_0 = DLDP_0_0 + DLDP_0_1*X_ref + DLDP_0_2*Y_ref
	DLDP_1 = DLDP_1_0 + DLDP_1_1*X_ref + DLDP_1_2*Y_ref + DLDP_1_3*X_ref**2 + DLDP_1_4*X_ref*Y_ref + DLDP_1_5*Y_ref**2

	return np.array([DLDP_0, DLDP_1]).transpose()

