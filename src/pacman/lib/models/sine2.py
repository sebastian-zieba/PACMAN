import sys

import numpy as np

sys.path.insert(0,'..')


def sine2(t, data, params, visit):
    a1, omega1, phi1, a2, omega2, phi2, a3, omega3, phi3, a12, omega12, phi12, a22, omega22, phi22, a32, omega32, phi32 = params

    a1 = a1[visit]
    omega1 = omega1[visit]
    phi1 = phi1[visit]
    a2 = a2[visit]
    omega2 = omega2[visit]
    phi2 = phi2[visit]
    a3 = a3[visit]
    omega3 = omega3[visit]
    phi3 = phi3[visit]

    a12 = a12[visit]
    omega12 = omega12[visit]
    phi12 = phi12[visit]
    a22 = a22[visit]
    omega22 = omega22[visit]
    phi22 = phi22[visit]
    a32 = a32[visit]
    omega32 = omega32[visit]
    phi32 = phi32[visit]

    #print "a1, omega1, phi1, a2, omega2, a3, omega3, phi3, phi2, a12, omega12, phi12, a22, omega22, phi22, a32, omega32, phi32",  a1, omega1, phi1, a2, omega2, a3, omega3, phi3, phi2, a12, omega12, phi12, a22, omega22, phi22, a32, omega32, phi32 

    TWOPI = np.pi*2.

    t_vis = data.t_vis[data.vis_idx[visit]]

    if t[0] < 2456300.:
        return ( 1. + a1*np.sin(TWOPI*omega1*t_vis) + phi1*np.cos(TWOPI*omega1*t_vis)
                    + a2*np.sin(TWOPI*omega2*t_vis) + phi2*np.cos(TWOPI*omega2*t_vis) 
                    + a3*np.sin(TWOPI*omega3*t_vis) + phi3*np.cos(TWOPI*omega3*t_vis)
               )
    else:
        return ( 1. + a12*np.sin(TWOPI*omega12*t_vis) + phi12*np.cos(TWOPI*omega12*t_vis) 
                    + a22*np.sin(TWOPI*omega22*t_vis) + phi22*np.cos(TWOPI*omega22*t_vis) 
                    + a32*np.sin(TWOPI*omega32*t_vis) + phi32*np.cos(TWOPI*omega32*t_vis) 
               )

           
