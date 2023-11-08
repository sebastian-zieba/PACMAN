            elif f == "polynomial2":
                self.sys.append(polynomial2)
                self.sys_porder.append([
                    data.par_order['v']*data.nvisit,
                    data.par_order['v2']*data.nvisit
                ])
            elif f == "polynomial3":
                self.sys.append(polynomial3)
                self.sys_porder.append([
                    data.par_order['v']*data.nvisit,
                    data.par_order['v2']*data.nvisit,
                    data.par_order['v3']*data.nvisit
                ])
            elif f == "logarithmic_visit":
                self.sys.append(logarithmic_visit)
                self.sys_porder.append([
                    data.par_order['log1']*data.nvisit,
                    data.par_order['log2']*data.nvisit
                ])
