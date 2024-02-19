import batman
import numpy as np


def transit(t, data, params, visit = 0):
    p = batman.TransitParams()

    t0, per, rp, a, inc, ecc, w, u1, u2 = params
    if data.ld_model == 2: p.limb_dark = "quadratic"
    elif data.ld_model == 1: p.limb_dark = "linear"
    elif data.ld_model == 'kipping2013': p.limb_dark = "quadratic"
    else: 
        print("unsupported limb darkening parameter")
        return 0
    #print(params)
    p.t0 = t0[visit] + data.toffset
    p.per = per[visit]
    p.rp = rp[visit]
    p.a = a[visit]
    p.inc = inc[visit]
    p.ecc = ecc[visit]
    p.w = w[visit]
    if data.ld_model == 2:
        p.u = [u1[visit], u2[visit]]
    elif data.ld_model == 1:
        p.u = [u1[visit]]
    elif data.ld_model == 'kipping2013':
        u1_quad = 2 * np.sqrt(u1[visit]) * u2[visit] # Eq. 15 in Kipping 2013
        u2_quad = np.sqrt(u1[visit]) * (1 - 2 * u2[visit])  # Eq. 16 in Kipping 2013
        p.u = np.array([u1_quad, u2_quad])
    m = batman.TransitModel(
        p, t, supersample_factor=5, exp_time = data.exp_time/24./60./60.
    )
    return m.light_curve(p)
