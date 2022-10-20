import batman


def transit(t, data, params, visit = 0):
    p = batman.TransitParams()

    t0, per, rp, a, inc, ecc, w, u1, u2 = params
    if data.ld_model == 2: p.limb_dark = "quadratic"
    elif data.ld_model == 1: p.limb_dark = "linear"
    else: 
        print("unsupported limb darkening parameter")
        return 0

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
    m = batman.TransitModel(
        p, t, supersample_factor=3, exp_time = data.exp_time/24./60./60.
    )
    return m.light_curve(p)
