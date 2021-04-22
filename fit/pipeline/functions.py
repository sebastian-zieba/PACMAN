import sys
sys.path.insert(0, './models')
from constant import constant
from polynomial1 import polynomial1
from polynomial2 import polynomial2
from sine1 import sine1
from sine2 import sine2
from upstream_downstream import upstream_downstream
from transit import transit
from eclipse import eclipse
from model_ramp import model_ramp
from divide_white import divide_white
from ackbar import ackbar

#need to automate appending parameters to functions
class Functions:
    def __init__(self, data, funcs):
        self.astro = []
        self.astro_porder = []
        self.sys = []
        self.sys_porder = []

        for f in funcs:
            if f == "constant":
                self.sys.append(constant)
                self.sys_porder.append([data.par_order['c']*data.nvisit]) #:(1 + data.par_order['c'])*data.nvisit])
            elif f == "upstream_downstream":
                self.sys.append(upstream_downstream)
                self.sys_porder.append([data.par_order['scale']*data.nvisit])
            elif f == "polynomial1":
                self.sys.append(polynomial1)
                self.sys_porder.append([data.par_order['v']*data.nvisit])
            elif f == "polynomial2":
                self.sys.append(polynomial2)
                self.sys_porder.append([
                    data.par_order['v']*data.nvisit,
                    data.par_order['v2']*data.nvisit
                ]) 
            elif f == "sine1":
                self.sys.append(sine1)
                self.sys_porder.append([
                    data.par_order['a1']*data.nvisit,
                    data.par_order['omega1']*data.nvisit,
                    data.par_order['phi1']*data.nvisit
                ]) 
            elif f == "sine2":
                self.sys.append(sine2)
                self.sys_porder.append([
                    data.par_order['a1']*data.nvisit,
                    data.par_order['omega1']*data.nvisit,
                    data.par_order['phi1']*data.nvisit,
                    data.par_order['a2']*data.nvisit,
                    data.par_order['omega2']*data.nvisit,
                    data.par_order['phi2']*data.nvisit,
                    data.par_order['a3']*data.nvisit,
                    data.par_order['omega3']*data.nvisit,
                    data.par_order['phi3']*data.nvisit,
                    data.par_order['a12']*data.nvisit,
                    data.par_order['omega12']*data.nvisit,
                    data.par_order['phi12']*data.nvisit,
                    data.par_order['a22']*data.nvisit,
                    data.par_order['omega22']*data.nvisit,
                    data.par_order['phi22']*data.nvisit,
                    data.par_order['a32']*data.nvisit,
                    data.par_order['omega32']*data.nvisit,
                    data.par_order['phi32']*data.nvisit
                ]) 
            elif f == "model_ramp":
                self.sys.append(model_ramp)
                self.sys_porder.append([
                    data.par_order['r1']*data.nvisit,
                    data.par_order['r2']*data.nvisit,
                    data.par_order['r3']*data.nvisit
                ]) 
            elif f == "ackbar":
                self.sys.append(ackbar)
                self.sys_porder.append([
                    data.par_order['trap_pop_s']*data.nvisit,
                    data.par_order['trap_pop_f']*data.nvisit,
                    data.par_order['dTrap_s']*data.nvisit,
                    data.par_order['dTrap_f']*data.nvisit
                ]) 
            elif f == "divide_white":
                self.sys.append(divide_white)
                self.sys_porder.append([])
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
                    data.par_order['u2']*data.nvisit, 
                    data.par_order['limb_dark']*data.nvisit
                    ])
            elif f == "eclipse":
                self.astro.append(eclipse)
                self.astro_porder.append([
                    data.par_order['t_secondary']*data.nvisit,
                    data.par_order['per']*data.nvisit, 
                    data.par_order['rp']*data.nvisit, 
                    data.par_order['fp']*data.nvisit, 
                    data.par_order['a']*data.nvisit, 
                    data.par_order['inc']*data.nvisit, 
                    data.par_order['ecc']*data.nvisit, 
                    data.par_order['w']*data.nvisit, 
                    ])
            else:
                #FIXME return error here
                return 0


    #def modelramp(self, t, params):
                    
