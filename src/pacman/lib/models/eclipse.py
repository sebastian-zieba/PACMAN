import batman


def eclipse(t, data, params, visit: float = 0.):
    p = batman.TransitParams()

    t_secondary, per, rp, fp, a, inc, ecc, w = params
    t_secondary = t_secondary[visit]
    per = per[visit]
    rp = rp[visit]
    fp = fp[visit]
    a = a[visit]
    inc = inc[visit]
    ecc = ecc[visit]
    w = w[visit]

    p.t_secondary = t_secondary
    p.per = per
    p.rp = rp
    p.fp = fp
    p.a = a
    p.inc = inc
    p.ecc = ecc
    p.w = w
    #TODO: WHY IS THIS HARD CODED?
    p.limb_dark = 'quadratic'
    p.u = [0.1, 0.2]
    
    m = batman.TransitModel(
        p, t, transittype = "secondary", supersample_factor=3, exp_time = data.exp_time/24./60./60.
    )
    return m.light_curve(p)
