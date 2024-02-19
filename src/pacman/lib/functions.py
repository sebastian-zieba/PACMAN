import sys
sys.path.insert(0, './models')
from ..lib.models.constant import constant
from ..lib.models.polynomial1 import polynomial1
from ..lib.models.polynomial1_full import polynomial1_full
from ..lib.models.polynomial2 import polynomial2
from ..lib.models.polynomial2_full import polynomial2_full
from ..lib.models.logarithmic_visit import logarithmic_visit
from ..lib.models.logarithmic_visit_full import logarithmic_visit_full
from ..lib.models.exponential_visit import exponential_visit
from ..lib.models.exponential_visit_full import exponential_visit_full
from ..lib.models.sine1 import sine1
from ..lib.models.sine2 import sine2
from ..lib.models.upstream_downstream import upstream_downstream
from ..lib.models.transit import transit
from ..lib.models.eclipse import eclipse
from ..lib.models.model_ramp import model_ramp
from ..lib.models.divide_white import divide_white
from ..lib.models.ackbar import ackbar
from ..lib.models.gp_sho import gp_sho
from ..lib.models.gp_matern32 import gp_matern32
from ..lib.models.constants_cj import constants_cj
from ..lib.models.uncmulti import uncmulti


#need to automate appending parameters to functions
class Functions:
    def __init__(self, data, funcs):
        self.astro = []
        self.astro_porder = []
        self.sys = []
        self.sys_porder = []
        self.gp = []
        self.gp_porder = []

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
            elif f == "polynomial1_full":
                self.sys.append(polynomial1_full)
                self.sys_porder.append([
                    data.par_order['scale']*data.nvisit,
                    data.par_order['c']*data.nvisit,
                    data.par_order['v']*data.nvisit
                ])
            elif f == "polynomial2":
                self.sys.append(polynomial2)
                self.sys_porder.append([
                    data.par_order['v']*data.nvisit,
                    data.par_order['v2']*data.nvisit
                ])
            elif f == "polynomial2_full":
                self.sys.append(polynomial2_full)
                self.sys_porder.append([
                    data.par_order['scale']*data.nvisit,
                    data.par_order['c']*data.nvisit,
                    data.par_order['v']*data.nvisit,
                    data.par_order['v2']*data.nvisit
                ])
            elif f == "logarithmic_visit":
                self.sys.append(logarithmic_visit)
                self.sys_porder.append([
                    data.par_order['log1']*data.nvisit,
                    data.par_order['log2']*data.nvisit
                ])
            elif f == "logarithmic_visit_full":
                self.sys.append(logarithmic_visit_full)
                self.sys_porder.append([
                    data.par_order['scale']*data.nvisit,
                    data.par_order['c']*data.nvisit,
                    data.par_order['log1']*data.nvisit,
                    data.par_order['log2']*data.nvisit
                ])
            elif f == "exponential_visit":
                self.sys.append(exponential_visit)
                self.sys_porder.append([
                    data.par_order['exp1']*data.nvisit,
                    data.par_order['exp2']*data.nvisit
                ])
            elif f == "exponential_visit_full":
                self.sys.append(exponential_visit_full)
                self.sys_porder.append([
                    data.par_order['scale']*data.nvisit,
                    data.par_order['c']*data.nvisit,
                    data.par_order['exp1']*data.nvisit,
                    data.par_order['exp2']*data.nvisit
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
                    data.par_order['u2']*data.nvisit
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
                    data.par_order['w']*data.nvisit
                    ])
            elif f == "gp_sho":
                self.gp.append(gp_sho)
                self.gp_porder.append([
                    data.par_order['logQ_gp']*data.nvisit,
                    data.par_order['logw_gp']*data.nvisit,
                    data.par_order['logS_gp']*data.nvisit,
                    data.par_order['log_jit']*data.nvisit
                    ])
            elif f == "gp_matern32":
                self.gp.append(gp_matern32)
                self.gp_porder.append([
                    data.par_order['log_rho']*data.nvisit,
                    data.par_order['log_sigma']*data.nvisit
                    ])
            elif f == "constants_cj":
                self.sys.append(constants_cj)
                cj_list = []
                for i in range(data.imax):
                    ci = f'c{i}'
                    cj_list.append(data.par_order[ci] * data.nvisit)
                self.sys_porder.append(cj_list)
            elif f == "uncmulti":
                self.sys.append(uncmulti)
                self.sys_porder.append([
                    data.par_order['uncmulti_val'] * data.nvisit
                ])
            else:
                #FIXME return error here
                return 0


    #def modelramp(self, t, params):
                    
