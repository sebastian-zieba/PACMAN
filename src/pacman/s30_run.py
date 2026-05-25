import time
import shutil
from pathlib import Path

import numpy as np
from astropy.io import ascii

from .lib import manageevent as me
from .lib import nice_fit_par
from .lib import plots
from .lib import util
from .lib.formatter import ReturnParams
from .lib.least_squares import lsq_fit
from .lib.mcmc import mcmc_fit
from .lib.model import Model
from .lib.nested import nested_sample
from .lib.read_data import Data
from .lib import logedit
from .lib import read_pcf as rd
from .lib import read_fit_par


def run30(pcf_path: Path, meta=None):
    """
    This functions reads in the spectroscopic or white light curve(s) and
    fits a model to them.
    """
    pcf_path = Path(pcf_path)

    # Read live/current obs_par.pcf to decide whether Stage 30 uses
    # Stage 20 white light curves or Stage 21 spectroscopic light curves.
    pcf = rd.read_pcf(pcf_path / "obs_par.pcf")
    fit_white = pcf.s30_fit_white.get(0)
    fit_spec = pcf.s30_fit_spec.get(0)

    if fit_white and fit_spec:
        raise ValueError("Only one of s30_fit_white or s30_fit_spec can be True.")

    if not fit_white and not fit_spec:
        raise ValueError("Either s30_fit_white or s30_fit_spec must be True.")

    if fit_spec:
        previous_stage_num = "21"
        stage_subdir = "spec_lc"
        copy_extracted_sp = True
        copy_extracted_lc = False
    else:
        previous_stage_num = "20"
        stage_subdir = "white_lc"
        copy_extracted_sp = False
        copy_extracted_lc = True

    meta, log = util.setup_stage(
        pcf_path=pcf_path,
        stage_num="30",
        previous_stage_num=previous_stage_num,
        stage_subdir=stage_subdir,
        copy_filelist=True,
        copy_xrefyref=True,
        copy_ancil=True,
        copy_extracted_lc=copy_extracted_lc,
        copy_extracted_sp=copy_extracted_sp,
        meta=meta,
    )

    # Create a directory for the white or the spectroscopic fit.
    # The Stage 30 run directory is already timestamped, so no extra
    # timestamped fit_* subdirectory is needed.
    if meta.s30_fit_white:
        meta.fitdir = Path("fit_white")
    elif meta.s30_fit_spec:
        meta.fitdir = Path("fit_spec")

    fit_dir = meta.workdir / meta.fitdir
    fit_dir.mkdir(parents=True, exist_ok=True)

    # Make fit_par nicer by lining up the columns
    nice_fit_par.nice_fit_par(meta.workdir / "fit_par.txt")

    # Reads in fit parameters from the fit_par file
    fit_par = ascii.read(
        meta.workdir / "fit_par.txt",
        format="commented_header",
        delimiter=r"\s",
        guess=False,
        fast_reader=False,
        fill_values=[("", "0")],
    )
    read_fit_par.validate_fit_par(
        fit_par,
        run_mcmc=meta.run_mcmc,
        run_nested=meta.run_nested,
        nsigma_lsq=5.0,
        warn_only_for_x=False, # I set that to False because for nested sampling the code will definitely break when "X" is used as a prior. 
    )

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

                read_fit_par.validate_c_against_light_curve(
                    fit_par,
                    data,
                    nsigma_lc=100.0,
                    nsigma_prior=100.0,
                    warn_only=False,
                )

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

                    read_fit_par.validate_c_against_light_curve(
                        fit_par,
                        data,
                        nsigma_lc=100.0,
                        nsigma_prior=100.0,
                        warn_only=False,
                    )

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

            if meta.run_verbose:
                print("rms, chi2red = ", model.rms, model.chi2red)

            # Save white systematics file if it was a white fit
            if meta.s30_fit_white:
                with (fit_dir / 'white_systematics.txt')\
                        .open("w", encoding="ascii") as outfile:
                    for i in range(len(model.all_sys)):
                        print(model.all_sys[i], file=outfile)
                    print('Saved white_systematics.txt file')

        meta.labels = data.free_parnames

        if meta.run_mcmc:
            print('\n*STARTS MCMC*')
            time.sleep(1.01)
            if meta.rescale_uncert:
                # NOTE: Rescale error bars so reduced chi-squared is one
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

    log.writelog("Saving Metadata")
    me.saveevent(meta, meta.workdir / "WFC3_Meta_Save", save=[])

    log.writelog("Finished s30")
    log.closelog()
    return meta


def update_clips(clips_array):
    clips_old = clips_array[0]
    clips_new = clips_array[1]
    clips_new[0] + sum([i <= clips_new[0] for i in clips_old])
    clips_new_updated = [clips_new[ii] + sum([i <= clips_new[ii] for i in clips_old]) for ii in range(len(clips_new))]
    return np.concatenate((clips_old, clips_new_updated))
