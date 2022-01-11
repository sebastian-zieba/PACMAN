import numpy as np


class FormatParams:
    """
    doc
    """
    def __init__(self, params, data):
        self.per = params[data.par_order['per']*data.nvisit:(1 + data.par_order['per'])*data.nvisit]
        self.t0 = params[data.par_order['t0']*data.nvisit:(1 + data.par_order['t0'])*data.nvisit]
        self.t_secondary = params[data.par_order['t_secondary']*data.nvisit:(1 + data.par_order['t_secondary'])*data.nvisit]
        self.w = params[data.par_order['w']*data.nvisit:(1 + data.par_order['w'])*data.nvisit]
        self.a = params[data.par_order['a']*data.nvisit:(1 + data.par_order['a'])*data.nvisit]
        self.inc = params[data.par_order['inc']*data.nvisit:(1 + data.par_order['inc'])*data.nvisit]
        self.rp = params[data.par_order['rp']*data.nvisit:(1 + data.par_order['rp'])*data.nvisit]
        self.fp = params[data.par_order['fp']*data.nvisit:(1 + data.par_order['fp'])*data.nvisit]
        self.u1 = params[data.par_order['u1']*data.nvisit:(1 + data.par_order['u1'])*data.nvisit]
        self.u2 = params[data.par_order['u2']*data.nvisit:(1 + data.par_order['u2'])*data.nvisit]
        self.ecc = params[data.par_order['ecc']*data.nvisit:(1 + data.par_order['ecc'])*data.nvisit]
        self.limb_dark = params[data.par_order['limb_dark']*data.nvisit:(1 + data.par_order['limb_dark'])*data.nvisit]
        self.c = params[data.par_order['c']*data.nvisit:(1 + data.par_order['c'])*data.nvisit]
        self.v = params[data.par_order['v']*data.nvisit:(1 + data.par_order['v'])*data.nvisit]
        self.v2 = params[data.par_order['v2']*data.nvisit:(1 + data.par_order['v2'])*data.nvisit]
        self.r1 = params[data.par_order['r1']*data.nvisit:(1 + data.par_order['r1'])*data.nvisit]
        self.r2 = params[data.par_order['r2']*data.nvisit:(1 + data.par_order['r2'])*data.nvisit]
        self.r3 = params[data.par_order['r3']*data.nvisit:(1 + data.par_order['r3'])*data.nvisit]
        self.scale = params[data.par_order['scale']*data.nvisit:(1 + data.par_order['scale'])*data.nvisit]
        self.trap_pop_s = params[data.par_order['trap_pop_s']*data.nvisit:(1 + data.par_order['trap_pop_s'])*data.nvisit]
        self.trap_pop_f = params[data.par_order['trap_pop_f']*data.nvisit:(1 + data.par_order['trap_pop_f'])*data.nvisit]
        self.dTrap_s = params[data.par_order['dTrap_s']*data.nvisit:(1 + data.par_order['dTrap_s'])*data.nvisit]
        self.dTrap_f = params[data.par_order['dTrap_f']*data.nvisit:(1 + data.par_order['dTrap_f'])*data.nvisit]

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