import time
import shutil
from pathlib import Path

import numpy as np
from astropy.io import ascii

from .lib import manageevent as me
from .lib.read_data import Data
from .lib.model import Model
from .lib.least_squares import lsq_fit
from .lib.mcmc import mcmc_fit
from .lib.nested import nested_sample
from .lib.formatter import ReturnParams
from .lib import nice_fit_par
from .lib import plots
from .lib import util


def run30(eventlabel: str, workdir: Path, meta=None):
    """This functions reads in the spectroscopic or white light curve(s) and
    fits a model to them."""
    print('Starting s30')

    if meta is None:
        meta = me.loadevent(workdir / f'WFC3_{eventlabel}_Meta_Save')

    # Create directories for Stage 3 processing
    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')

    # Create a directory for the white or the spectroscopic fit
    if meta.s30_fit_white:
        meta.fitdir = Path('fit_white') / f'fit_{datetime}_{meta.eventlabel}'
    elif meta.s30_fit_spec:
        meta.fitdir = Path('fit_spec') / f'fit_{datetime}_{meta.eventlabel}'

    fit_dir = meta.workdir / meta.fitdir
    if not fit_dir.exists():
        fit_dir.mkdir(parents=True, exist_ok=True)

    # Make fit_par nicer by lining up the columns
    nice_fit_par.nice_fit_par(meta.workdir / "fit_par.txt")

    # Copy pcf and fit_par files into the new fitdirectory
    shutil.copy(meta.workdir / "obs_par.pcf", fit_dir)
    shutil.copy(meta.workdir / "fit_par.txt", fit_dir)

    # Reads in fit parameters from the fit_par file
    #TODO: Check that fit_par is configured correctly. Eg initial value has to be within boundaries!
    fit_par = ascii.read(meta.workdir / "fit_par.txt", Reader=ascii.CommentedHeader)

    # Read in the user wanted fit functions
    myfuncs = meta.s30_myfuncs

    # Indentify path of light curve files
    files, meta = util.read_fitfiles(meta)

    # Prepare some lists for diagnostics
    # TODO: Maybe make dict out of this?
    if meta.run_verbose:
        vals = []
        errs = []
        idxs = []
    if meta.run_mcmc:
        vals_mcmc       = []
        errs_lower_mcmc = []
        errs_upper_mcmc = []
    if meta.run_nested:
        vals_nested = []
        errs_lower_nested = []
        errs_upper_nested = []

    meta = util.log_run_setup(meta)

    # Loop through light curves
    for counter, f in enumerate(files):
        print('\n****** File: {0}/{1}'.format(counter+1, len(files)))
        meta.s30_file_counter = counter
        time.sleep(1.01) #sleep to prevent overwriting of data if saved in the same second
        meta.run_file = f
        meta.fittime = time.strftime('%Y-%m-%d_%H-%M-%S')

        ### RUN LEAST SQUARES ###
        ## IF NO CLIPPING IS WISHED
        if meta.run_lsq:
            if meta.run_clipiters == 0:
                print('\n')
                data = Data(f, meta, fit_par)
                model = Model(data, myfuncs)
                print('\n*STARTS LEAST SQUARED*')
                data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True) #not clipping

            ## IF THE USER WANTS CLIPPING
            else:
                clip_idxs = []
                for iii in range(meta.run_clipiters+1):
                    print('\n')
                    print('Sigma Iters: ', iii, 'of', meta.run_clipiters)
                    if iii == 0:
                        data = Data(f, meta, fit_par)
                    else:
                        data = Data(f, meta, fit_par, clip_idx)
                    model = Model(data, myfuncs)
                    if iii == meta.run_clipiters:
                        data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
                    else:
                        data, model, params, clip_idx, m = lsq_fit(fit_par, data, meta, model, myfuncs)
                        print("rms, chi2red = ", model.rms, model.chi2red)
                        print(clip_idx == [])
                        if clip_idx == []: break
                        clip_idxs.append(clip_idx)
                        print(clip_idxs)
                        print('length: ', len(clip_idxs) )
                        if len(clip_idxs)>1:
                            clip_idx = update_clips(clip_idxs)
                            print(clip_idx)
                            clip_idxs = update_clips(clip_idxs)
                            print(clip_idxs)

            if meta.run_verbose == True:
                print("rms, chi2red = ", model.rms, model.chi2red)

            # Save white systematics file if it was a white fit
            if meta.s30_fit_white == True:
                outfile = open(fit_dir / 'white_systematics.txt', "w")
                for i in range(len(model.all_sys)): print(model.all_sys[i], file = outfile)
                print('Saved white_systematics.txt file')
                outfile.close()

        meta.labels = data.free_parnames

        if meta.run_mcmc:
            print('\n*STARTS MCMC*')
            time.sleep(1.01)
            if meta.rescale_uncert:
                ##rescale error bars so reduced chi-squared is one
                print(f'rescale_uncert in the pcf was set to {meta.rescale_uncert}')
                if model.chi2red < 1:
                    print('After the first fit, you got chi2_red < 1')
                    print(f'rescale_uncert_smaller_1 was set to {meta.rescale_uncert_smaller_1}')
                    if meta.rescale_uncert_smaller_1:
                        # also rescale if chi2red < 1
                        print('We will therefore rescale the errorbars, so that chi2_red = 1')
                        data.err *= np.sqrt(model.chi2red)
                    else:
                        print('chi2red < 1 and rescale_uncert_smaller_1 = False --> no rescaling')
                elif model.chi2red >= 1:
                    print('After the first fit, you got chi2_red >= 1')
                    print('Errorbars are being rescaled so that chi2_red = 1')
                    data.err *= np.sqrt(model.chi2red)
            # if a least square hasnt been run yet or the uncertainies should be rescaled, to the lsq again
            if not meta.run_lsq or meta.rescale_uncert:
                data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
                if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
            val_mcmc, err_lower_mcmc, err_upper_mcmc, fit = mcmc_fit(data, model, params, f, meta, fit_par)

        if meta.run_nested:
            print('\n*STARTS NESTED SAMPLING*')
            time.sleep(1.01)
            if meta.rescale_uncert:
                ##rescale error bars so reduced chi-squared is one
                print('Errorbars are being rescaled so that chi2_red = 1')
                if model.chi2red > 1:
                    data.err *= np.sqrt(model.chi2red)
            if not meta.run_lsq or meta.rescale_uncert:
                data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
                if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
            val_nested, err_lower_nested, err_upper_nested, fit = nested_sample(data, model, params, f, meta, fit_par)

        if meta.run_verbose:
            val, err, idx = ReturnParams(m, data)
            vals.append(val)
            errs.append(err)
            idxs.append(idx)

        if meta.run_mcmc:
            vals_mcmc.append(val_mcmc)
            errs_lower_mcmc.append(err_lower_mcmc)
            errs_upper_mcmc.append(err_upper_mcmc)

        if meta.run_nested:
            vals_nested.append(val_nested)
            errs_lower_nested.append(err_lower_nested)
            errs_upper_nested.append(err_upper_nested)

    if meta.run_verbose:
        plots.params_vs_wvl(vals, errs, idxs, meta)
        if meta.run_mcmc:
            plots.params_vs_wvl_mcmc(vals_mcmc, errs_lower_mcmc, errs_upper_mcmc, meta)

        if meta.run_nested:
            plots.params_vs_wvl_nested(vals_nested, errs_lower_nested, errs_upper_nested, meta)

        if not meta.s30_fit_white and ('rp' in meta.labels):
            # Saves rprs and wvl as a txt file
            util.make_lsq_rprs_txt(vals, errs, idxs, meta)
            # Saves rprs vs wvl as a plot
            plots.lsq_rprs(vals, errs, idxs, meta)

        if not meta.s30_fit_white and ('rp' in meta.labels) and meta.run_mcmc:
            plots.mcmc_rprs(vals_mcmc, errs_lower_mcmc, errs_upper_mcmc, meta)
            util.make_rprs_txt(vals_mcmc, errs_lower_mcmc, errs_upper_mcmc, meta, fitter='mcmc')

        if not meta.s30_fit_white and ('rp' in meta.labels) and meta.run_nested:
            plots.nested_rprs(vals_nested, errs_lower_nested, errs_upper_nested, meta)
            util.make_rprs_txt(vals_nested, errs_lower_nested, errs_upper_nested, meta, fitter='nested')

    if meta.run_mcmc or meta.run_nested:
        util.save_fit_output(fit, data, meta)

    print('Finished s30')

    return meta


def update_clips(clips_array):
    clips_old = clips_array[0]
    clips_new = clips_array[1]
    clips_new[0] + sum([i <= clips_new[0] for i in clips_old])
    clips_new_updated = [clips_new[ii] + sum([i <= clips_new[ii] for i in clips_old]) for ii in range(len(clips_new))]
    return np.concatenate((clips_old, clips_new_updated))
