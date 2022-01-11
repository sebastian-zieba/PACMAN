import numpy as np

#calculates the slope and intercept for the trace, given the position  of the direct image in physical pixels
#these coefficients are for the WFC3 G141 grism
def trace(X_ref, Y_ref):		
	BEAMA_i = 41 
	BEAMA_f = 248

	DYDX_0_0 = -3.55018E-01 
	DYDX_0_1 = 3.28722E-05
	DYDX_0_2 = -1.44571E-03

	DYDX_1_0 = 1.42852E-02 
	DYDX_1_1 = -7.20713E-06
	DYDX_1_2 = -2.42542E-06
	DYDX_1_3 = 1.18294E-09
	DYDX_1_4 = 1.19634E-08
	DYDX_1_5 = 6.17274E-10


	DYDX_0 = DYDX_0_0 + DYDX_0_1*X_ref + DYDX_0_2*Y_ref
	DYDX_1 = DYDX_1_0 + DYDX_1_1*X_ref + DYDX_1_2*Y_ref + DYDX_1_3*X_ref**2 + DYDX_1_4*X_ref*Y_ref + DYDX_1_5*Y_ref**2
	
	return [DYDX_0, DYDX_1]

#calculates coefficients for the dispersion solution
def dispersion(X_ref, Y_ref):		#X_ref and Y_ref are the centroid position in physical pixels
	DLDP_0_0 = 6.38738E+03 
	DLDP_0_1 = 4.55507E-02
	DLDP_0_2 = 0.0

	DLDP_1_0 = 2.35716E+01 
	DLDP_1_1 = 3.60396E-04
	DLDP_1_2 = 1.58739E-03
	DLDP_1_3 = -4.25234E-07
	DLDP_1_4 = -6.53726E-08
	DLDP_1_5 = -6.75872E-08

	DLDP_0 = DLDP_0_0 + DLDP_0_1*X_ref + DLDP_0_2*Y_ref
	DLDP_1 = DLDP_1_0 + DLDP_1_1*X_ref + DLDP_1_2*Y_ref + DLDP_1_3*X_ref**2 + DLDP_1_4*X_ref*Y_ref + DLDP_1_5*Y_ref**2

	return np.array([DLDP_0, DLDP_1]).transpose()

