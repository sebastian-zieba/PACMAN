        elif 'polynomial2' in data.s30_myfuncs:
            self.v = params[data.par_order['v']*data.nvisit:(1 + data.par_order['v'])*data.nvisit]
            self.v2 = params[data.par_order['v2']*data.nvisit:(1 + data.par_order['v2'])*data.nvisit]
        elif 'polynomial3' in data.s30_myfuncs:
            self.v = params[data.par_order['v']*data.nvisit:(1 + data.par_order['v'])*data.nvisit]
            self.v2 = params[data.par_order['v2']*data.nvisit:(1 + data.par_order['v2'])*data.nvisit]
            self.v3 = params[data.par_order['v3']*data.nvisit:(1 + data.par_order['v3'])*data.nvisit]
        elif 'logarithmic_visit' in data.s30_myfuncs:
            self.log1 = params[data.par_order['log1']*data.nvisit:(1 + data.par_order['log1'])*data.nvisit]
            self.log2 = params[data.par_order['log2']*data.nvisit:(1 + data.par_order['log2'])*data.nvisit]
