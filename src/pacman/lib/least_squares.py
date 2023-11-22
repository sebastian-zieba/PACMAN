import numpy as np
from astropy.stats import sigma_clip

from . import mpfit
from . import plots
from . import util
from . import read_fit_par
from .formatter import PrintParams
from .options import OPTIONS


def residuals(params, data, model, fjac=None):
    """Calculates the residuals of the fit."""
    fit = model.fit(data, params)
    return [0, fit.resid/data.err]


def lsq_fit(fit_par, data, meta, model, myfuncs, noclip=False):
    """Runs the least square routine using MPFIT."""
    # NOTE: Creates the lsq directory
    util.create_res_dir(meta)

    # TODO: noclip = True should be standard
    nvisit = data.nvisit
    npar = len(data.parnames)*nvisit

    # initializes least squares fit parameters
    parinfo = [{'value': 0, 'fixed': 0, 'limited': [0, 0],
                'limits': [0.0, 0.0], 'step': 0.0} for j in range(npar)]
    params_s = []

    # loops through parameters and visits
    # sets initial guess, step size, tie, bounds
    parinfo, params_s = read_fit_par.read_fit_par_for_ls(parinfo, params_s, data, fit_par)

    if meta.save_raw_lc_plot:
        plots.plot_raw(data, meta)
    fa = {'data': data, 'model': model}

    if 'divide_white' in data.s30_myfuncs and meta.s30_fit_spec:
        sys_vector = np.genfromtxt(meta.white_sys_path)
        data.all_sys = sys_vector
        # print("subtracting 2 from dof for divide-white")
        # data.nfree_param -= 2
        # data.dof += 2

    print('Runs MPFIT... ')
    m = mpfit.mpfit(residuals, params_s,
                    functkw=fa, parinfo=parinfo, quiet=True)

    if not noclip:
        # If user wants to sigma clip but there is nothing to clip:
        if sum(np.ma.getmask(sigma_clip(model.resid, sigma=meta.run_clipsigma, maxiters=1))) == 0:
            clip_idx = []
            if meta.save_fit_lc_plot: plots.plot_fit_lc2(data, model, meta)
        else:
            # If user wants to sigma clip and there are outliers:
            clip_idx = np.where(np.ma.getmask(sigma_clip(model.resid, sigma=meta.run_clipsigma, maxiters=1))==True)[0]
            print('Outlier Identified: ', len(clip_idx))
            print('Outlier idx: ', clip_idx)
            if meta.save_fit_lc_plot: plots.plot_fit_lc(data, model, meta)

    if m.errmsg: print("MPFIT error message", m.errmsg)

    if meta.run_verbose:
        with (meta.workdir / meta.fitdir / 'lsq_res' /
              f"lsq_res_bin{meta.s30_file_counter}_wvl{meta.wavelength:0.3f}.txt")\
                      .open("w", encoding=OPTIONS["encoding"]) as f_lsq:
            PrintParams(m, data, savefile=f_lsq)
        PrintParams(m, data)

    if meta.save_fit_lc_plot:
        if not (meta.workdir / meta.fitdir / 'fit_lc').exists():
            (meta.workdir / meta.fitdir / 'fit_lc').mkdir(parents=True)
        plots.plot_fit_lc2(data, model, meta)
        plots.plot_fit_lc3(data, model, meta)
        plots.save_plot_raw_data(data, meta)
        plots.save_astrolc_data(data, model, meta)

    util.append_fit_output(model, meta, fitter='lsq')

    if meta.save_allan_plot:
        plots.rmsplot(model, data, meta, fitter='lsq')

    if noclip:
        return data, model, m.params, m
    return data, model, m.params, clip_idx, m
