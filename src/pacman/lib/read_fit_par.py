import numpy as np


def read_fit_par_for_ls(parinfo, params_s, data, fit_par):
    """
    This function reads in the rows in the fit_par.txt file and saves it into a format so that it can be used with MPFIT.
    """
    nvisit = data.nvisit
    ii = 0
    for i in range(int(len(data.parnames))):
        if str(fit_par['tied'][ii]) == "-1":
            for j in range(nvisit):
                parinfo[i*nvisit+j]['value'] = fit_par['value'][ii]
                parinfo[i*nvisit+j]['step'] = fit_par['step_size'][ii]
                parinfo[i*nvisit+j]['fixed'] = fit_par['fixed'][ii].lower() == "true"
                if j>0 and str(fit_par['tied'][ii]) == "-1":
                    parinfo[i*nvisit+j]['tied'] = 'p[{0}]'.format(nvisit*i)
                if fit_par['lo_lim'][ii].lower() == "true":
                    parinfo[i*nvisit+j]['limited'][0] = True
                    parinfo[i*nvisit+j]['limits'][0] = fit_par['lo_val'][ii]
                if fit_par['hi_lim'][ii].lower() == "true":
                    parinfo[i*nvisit+j]['limited'][1] = True
                    parinfo[i*nvisit+j]['limits'][1] = fit_par['hi_val'][ii]
                params_s.append(fit_par['value'][ii])
            ii = ii + 1
        else:
            for j in range(nvisit):
                parinfo[i * nvisit + j]['value'] = fit_par['value'][ii]
                parinfo[i * nvisit + j]['step'] = fit_par['step_size'][ii]
                parinfo[i * nvisit + j]['fixed'] = fit_par['fixed'][ii].lower() == "true"
                if j > 0 and str(fit_par['tied'][ii]) == "-1":
                    parinfo[i * nvisit + j]['tied'] = 'p[{0}]'.format(nvisit * i)
                if fit_par['lo_lim'][ii].lower() == "true":
                    parinfo[i * nvisit + j]['limited'][0] = True
                    parinfo[i * nvisit + j]['limits'][0] = fit_par['lo_val'][ii]
                if fit_par['hi_lim'][ii].lower() == "true":
                    parinfo[i * nvisit + j]['limited'][1] = True
                    parinfo[i * nvisit + j]['limits'][1] = fit_par['hi_val'][ii]
                params_s.append(fit_par['value'][ii])
                ii = ii + 1
    return parinfo, np.array(params_s)


def get_step_size(data, params, meta, fit_par):
    """
    Get the step sizes which were set in the fit_par.txt file by the user
    """
    nvisit = int(meta.nvisit)
    nfree_param = data.nfree_param
    step_size = []
    ii = 0
    for i in range(len(fit_par)):
        if ii == len(fit_par):
            break
        if fit_par['fixed'][ii].lower() == "false":
                if str(fit_par['tied'][ii]) == "-1":
                    step_size.append(fit_par['step_size'][ii])
                    ii = ii + 1
                else:
                    for j in range(nvisit):
                        step_size.append(fit_par['step_size'][ii])
                        ii = ii + 1
        else:
            ii = ii + 1

    return np.array(step_size)