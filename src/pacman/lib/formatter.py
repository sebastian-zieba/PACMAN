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
        elif 'model_ramp' in data.s30_myfuncs:
            self.r1 = params[data.par_order['r1']*data.nvisit:(1 + data.par_order['r1'])*data.nvisit]
            self.r2 = params[data.par_order['r2']*data.nvisit:(1 + data.par_order['r2'])*data.nvisit]
            self.r3 = params[data.par_order['r3']*data.nvisit:(1 + data.par_order['r3'])*data.nvisit]
        elif 'upstream_downstream' in data.s30_myfuncs:
            self.scale = params[data.par_order['scale']*data.nvisit:(1 + data.par_order['scale'])*data.nvisit]
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
        elif 'constants_cj' in data.s30_myfuncs:
            self.c0 = params[data.par_order['c0']*data.nvisit:(1 + data.par_order['c0'])*data.nvisit]
            self.c1 = params[data.par_order['c1']*data.nvisit:(1 + data.par_order['c1'])*data.nvisit]
            self.c2 = params[data.par_order['c2']*data.nvisit:(1 + data.par_order['c2'])*data.nvisit]
            self.c3 = params[data.par_order['c3']*data.nvisit:(1 + data.par_order['c3'])*data.nvisit]
            self.c4 = params[data.par_order['c4']*data.nvisit:(1 + data.par_order['c4'])*data.nvisit]
            self.c5 = params[data.par_order['c5']*data.nvisit:(1 + data.par_order['c5'])*data.nvisit]
            self.c6 = params[data.par_order['c6']*data.nvisit:(1 + data.par_order['c6'])*data.nvisit]
            self.c7 = params[data.par_order['c7']*data.nvisit:(1 + data.par_order['c7'])*data.nvisit]
            self.c8 = params[data.par_order['c8']*data.nvisit:(1 + data.par_order['c8'])*data.nvisit]
            self.c9 = params[data.par_order['c9']*data.nvisit:(1 + data.par_order['c9'])*data.nvisit]
            self.c10 = params[data.par_order['c10']*data.nvisit:(1 + data.par_order['c10'])*data.nvisit]
            self.c11 = params[data.par_order['c11']*data.nvisit:(1 + data.par_order['c11'])*data.nvisit]
            self.c12 = params[data.par_order['c12']*data.nvisit:(1 + data.par_order['c12'])*data.nvisit]
            self.c13 = params[data.par_order['c13']*data.nvisit:(1 + data.par_order['c13'])*data.nvisit]
            self.c14 = params[data.par_order['c14']*data.nvisit:(1 + data.par_order['c14'])*data.nvisit]
            self.c15 = params[data.par_order['c15']*data.nvisit:(1 + data.par_order['c15'])*data.nvisit]
            self.c16 = params[data.par_order['c16']*data.nvisit:(1 + data.par_order['c16'])*data.nvisit]
            self.c17 = params[data.par_order['c17']*data.nvisit:(1 + data.par_order['c17'])*data.nvisit]
            self.c18 = params[data.par_order['c18']*data.nvisit:(1 + data.par_order['c18'])*data.nvisit]
            self.c19 = params[data.par_order['c19']*data.nvisit:(1 + data.par_order['c19'])*data.nvisit]
            self.c20 = params[data.par_order['c20']*data.nvisit:(1 + data.par_order['c20'])*data.nvisit]
            self.c21 = params[data.par_order['c21']*data.nvisit:(1 + data.par_order['c21'])*data.nvisit]
            self.c22 = params[data.par_order['c22']*data.nvisit:(1 + data.par_order['c22'])*data.nvisit]
            self.c23 = params[data.par_order['c23']*data.nvisit:(1 + data.par_order['c23'])*data.nvisit]
            self.c24 = params[data.par_order['c24']*data.nvisit:(1 + data.par_order['c24'])*data.nvisit]
            self.c25 = params[data.par_order['c25']*data.nvisit:(1 + data.par_order['c25'])*data.nvisit]
            self.c26 = params[data.par_order['c26']*data.nvisit:(1 + data.par_order['c26'])*data.nvisit]
            self.c27 = params[data.par_order['c27']*data.nvisit:(1 + data.par_order['c27'])*data.nvisit]
            self.c28 = params[data.par_order['c28']*data.nvisit:(1 + data.par_order['c28'])*data.nvisit]
            self.c29 = params[data.par_order['c29']*data.nvisit:(1 + data.par_order['c29'])*data.nvisit]


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