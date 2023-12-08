import numpy as np


class FormatParams:
    """
    doc
    """
    def __init__(self, params, data):
        if 'transit' in data.s30_myfuncs:
            self.per = params[data.par_order['per']*data.nvisit:(1 + data.par_order['per'])*data.nvisit]
            self.t0 = params[data.par_order['t0']*data.nvisit:(1 + data.par_order['t0'])*data.nvisit]
            self.w = params[data.par_order['w']*data.nvisit:(1 + data.par_order['w'])*data.nvisit]
            self.a = params[data.par_order['a']*data.nvisit:(1 + data.par_order['a'])*data.nvisit]
            self.inc = params[data.par_order['inc']*data.nvisit:(1 + data.par_order['inc'])*data.nvisit]
            self.rp = params[data.par_order['rp']*data.nvisit:(1 + data.par_order['rp'])*data.nvisit]
            self.u1 = params[data.par_order['u1'] * data.nvisit:(1 + data.par_order['u1']) * data.nvisit]
            self.u2 = params[data.par_order['u2'] * data.nvisit:(1 + data.par_order['u2']) * data.nvisit]
            self.ecc = params[data.par_order['ecc'] * data.nvisit:(1 + data.par_order['ecc']) * data.nvisit]
        elif 'eclipse' in data.s30_myfuncs:
            self.t_secondary = params[data.par_order['t_secondary']*data.nvisit:(1 + data.par_order['t_secondary'])*data.nvisit]
            self.w = params[data.par_order['w']*data.nvisit:(1 + data.par_order['w'])*data.nvisit]
            self.a = params[data.par_order['a']*data.nvisit:(1 + data.par_order['a'])*data.nvisit]
            self.inc = params[data.par_order['inc']*data.nvisit:(1 + data.par_order['inc'])*data.nvisit]
            self.rp = params[data.par_order['rp']*data.nvisit:(1 + data.par_order['rp'])*data.nvisit]
            self.u1 = params[data.par_order['u1'] * data.nvisit:(1 + data.par_order['u1']) * data.nvisit]
            self.u2 = params[data.par_order['u2'] * data.nvisit:(1 + data.par_order['u2']) * data.nvisit]
            self.ecc = params[data.par_order['ecc'] * data.nvisit:(1 + data.par_order['ecc']) * data.nvisit]
            self.fp = params[data.par_order['fp']*data.nvisit:(1 + data.par_order['fp'])*data.nvisit]
        elif 'constant' in data.s30_myfuncs:
            self.c = params[data.par_order['c']*data.nvisit:(1 + data.par_order['c'])*data.nvisit]
        elif 'polynomial1' in data.s30_myfuncs:
            self.v = params[data.par_order['v'] * data.nvisit:(1 + data.par_order['v']) * data.nvisit]
        elif 'polynomial2' in data.s30_myfuncs:
            self.v = params[data.par_order['v']*data.nvisit:(1 + data.par_order['v'])*data.nvisit]
            self.v2 = params[data.par_order['v2']*data.nvisit:(1 + data.par_order['v2'])*data.nvisit]
        elif 'logarithmic_visit' in data.s30_myfuncs:
            self.log1 = params[data.par_order['log1']*data.nvisit:(1 + data.par_order['log1'])*data.nvisit]
            self.log2 = params[data.par_order['log2']*data.nvisit:(1 + data.par_order['log2'])*data.nvisit]
        elif 'exponential_visit' in data.s30_myfuncs:
            self.exp1 = params[data.par_order['exp1']*data.nvisit:(1 + data.par_order['exp1'])*data.nvisit]
            self.exp2 = params[data.par_order['exp2']*data.nvisit:(1 + data.par_order['exp2'])*data.nvisit]
        elif 'model_ramp' in data.s30_myfuncs:
            self.r1 = params[data.par_order['r1']*data.nvisit:(1 + data.par_order['r1'])*data.nvisit]
            self.r2 = params[data.par_order['r2']*data.nvisit:(1 + data.par_order['r2'])*data.nvisit]
            self.r3 = params[data.par_order['r3']*data.nvisit:(1 + data.par_order['r3'])*data.nvisit]
        elif 'upstream_downstream' in data.s30_myfuncs:
            self.scale = params[data.par_order['scale']*data.nvisit:(1 + data.par_order['scale'])*data.nvisit]
        elif 'model_rowshift' in data.s30_myfuncs:
            self.rowshift_vf = params[data.par_order['rowshift_vf']*data.nvisit:(1 + data.par_order['rowshift_vf'])*data.nvisit]
            self.rowshift_vr = params[data.par_order['rowshift_vr']*data.nvisit:(1 + data.par_order['rowshift_vr'])*data.nvisit]
        elif 'ackbar' in data.s30_myfuncs:
            self.trap_pop_s = params[data.par_order['trap_pop_s']*data.nvisit:(1 + data.par_order['trap_pop_s'])*data.nvisit]
            self.trap_pop_f = params[data.par_order['trap_pop_f']*data.nvisit:(1 + data.par_order['trap_pop_f'])*data.nvisit]
            self.dTrap_s = params[data.par_order['dTrap_s']*data.nvisit:(1 + data.par_order['dTrap_s'])*data.nvisit]
            self.dTrap_f = params[data.par_order['dTrap_f']*data.nvisit:(1 + data.par_order['dTrap_f'])*data.nvisit]
        elif 'gp_sho' in data.s30_myfuncs:
            self.logS_gp = params[data.par_order['logS_gp'] * data.nvisit:(1 + data.par_order['logS_gp']) * data.nvisit]
            self.logw_gp = params[data.par_order['logw_gp'] * data.nvisit:(1 + data.par_order['logw_gp']) * data.nvisit]
            self.logQ_gp = params[data.par_order['logQ_gp'] * data.nvisit:(1 + data.par_order['logQ_gp']) * data.nvisit]
            self.log_jit = params[data.par_order['log_jit'] * data.nvisit:(1 + data.par_order['log_jit']) * data.nvisit]
        elif 'gp_matern32' in data.s30_myfuncs:
            self.log_rho = params[data.par_order['log_rho'] * data.nvisit:(1 + data.par_order['log_rho']) * data.nvisit]
            self.log_sigma = params[data.par_order['log_sigma'] * data.nvisit:(1 + data.par_order['log_sigma']) * data.nvisit]
        elif 'constants_cj' in data.s30_myfuncs:
            for i in range(data.imax):
                ci = f'c{i}'
                # TODO: So this might not be best practice (probably dicts would be better)
                #  But I think there's no way currently to make this better the way the code is written
                exec(f'self.c{i} = params[data.par_order[{ci}]*data.nvisit:(1 + data.par_order[{ci}])*data.nvisit]')
        elif 'uncmulti' in data.s30_myfuncs:
            self.uncmulti_val = params[data.par_order['uncmulti_val'] * data.nvisit:(1 + data.par_order['uncmulti_val']) * data.nvisit]


def PrintParams(m, data, savefile=False):
    for name in data.parnames:
        for vis in range(data.nvisit):
            if m.perror[data.par_order[name]*data.nvisit + vis] > 0.: 
                if not savefile:
                    print(name+"_"+str(vis), \
                          "\t", "{0:0.4e}".format(m.params[data.par_order[name]*data.nvisit + vis]), \
                          "\t", "{0:0.4e}".format(m.perror[data.par_order[name]*data.nvisit + vis]))
                else:
                    print(name+"_"+str(vis), \
                          "\t", "{0:0.4e}".format(m.params[data.par_order[name]*data.nvisit + vis]), \
                          "\t", "{0:0.4e}".format(m.perror[data.par_order[name]*data.nvisit + vis]), file=savefile)


def ReturnParams(m, data):
    val = []
    err = []
    for name in data.parnames:
        for vis in range(data.nvisit):
            if m.perror[data.par_order[name]*data.nvisit + vis] > 0.:
                val.append(m.params[data.par_order[name]*data.nvisit + vis])
                err.append(m.perror[data.par_order[name]*data.nvisit + vis])
            else:
                val.append(np.nan)
                err.append(np.nan)
    val = np.array(val)
    err = np.array(err)
    idx = np.arange(len(val), dtype=int)[~np.isnan(val)]
    return val, err, idx
