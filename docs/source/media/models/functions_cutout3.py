            elif f == "transit":
                self.astro.append(transit)
                self.astro_porder.append([
                    data.par_order['t0']*data.nvisit,
                    data.par_order['per']*data.nvisit, 
                    data.par_order['rp']*data.nvisit, 
                    data.par_order['a']*data.nvisit, 
                    data.par_order['inc']*data.nvisit, 
                    data.par_order['ecc']*data.nvisit, 
                    data.par_order['w']*data.nvisit, 
                    data.par_order['u1']*data.nvisit, 
                    data.par_order['u2']*data.nvisit
                    ])
