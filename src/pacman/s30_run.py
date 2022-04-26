import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from astropy.io import ascii
import getopt
import time
import shutil
import time as pythontime
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


def run30(eventlabel, workdir, meta=None):
    """
    This functions reads in the spectroscopic or white light curve(s) and fits a model to them.
    """
    print('Starting s30')

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    # Create directories for Stage 3 processing
    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
    meta.fitdir = '/fit_' + datetime + '_' + meta.eventlabel
    if not os.path.exists(meta.workdir + meta.fitdir):
        os.makedirs(meta.workdir + meta.fitdir)

    # Make fit_par nice...
    nice_fit_par.nice_fit_par(meta.workdir + "/fit_par.txt")

    # Copy pcf and fit_par files
    shutil.copy(meta.workdir + "/obs_par.pcf", meta.workdir + meta.fitdir)
    shutil.copy(meta.workdir + "/fit_par.txt", meta.workdir + meta.fitdir)

    # reads in fit parameters
    #TODO: Check that fit_par is configured correctly. Eg initial value has to be within boundaries!
    fit_par = ascii.read(meta.workdir + "/fit_par.txt", Reader=ascii.CommentedHeader)

    #read in the user wanted fit functions
    myfuncs = meta.s30_myfuncs

    #store files to fit
    files, meta = util.read_fitfiles(meta)

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

    meta.chi2red_list = []

    for counter, f in enumerate(files):
        print('\n****** File: {0}/{1}'.format(counter+1, len(files)))
        meta.s30_file_counter = counter
        time.sleep(1.1) #sleep to prevent overwriting of data if saved in the same second
        meta.run_file = f
        meta.fittime = time.strftime('%Y-%m-%d_%H-%M-%S')

        ### RUN LEAST SQUARES ###
        ## IF NO CLIPPING IS WISHED
        if meta.run_clipiters == 0:
            print('\n')
            data = Data(f, meta, fit_par)
            model = Model(data, myfuncs)
            print('*STARTS LEAST SQUARED*')
            data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True) #not clipping
            meta.chi2red_list.append(model.chi2red)
            if meta.save_fit_lc_plot: plots.plot_fit_lc2(data, model, meta)
            if meta.save_fit_lc_plot: plots.plot_fit_lc3(data, model, meta)
            if meta.save_fit_lc_plot: plots.save_plot_raw_data(data, meta)
            if meta.save_fit_lc_plot: plots.save_astrolc_data(data, model, meta)
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

        #FIXME : make this automatic!
        """outfile = open("white_systematics.txt", "w")
        for i in range(len(model.all_sys)): print(model.all_sys[i], file = outfile)
        outfile.close()"""

        meta.labels = labels_gen(params, meta, fit_par)

        if meta.run_mcmc:
            print('*STARTS MCMC*')
            time.sleep(1.1)
            if meta.rescale_uncert:
                ##rescale error bars so reduced chi-squared is one
                data.err *= np.sqrt(model.chi2red)
            data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
            if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
            val_mcmc, err_lower_mcmc, err_upper_mcmc = mcmc_fit(data, model, params, f, meta, fit_par)

        if meta.run_nested:
            time.sleep(1.1)
            if meta.rescale_uncert:
                ##rescale error bars so reduced chi-squared is one
                data.err *= np.sqrt(model.chi2red)
            data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
            if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
            val_nested, err_lower_nested, err_upper_nested = nested_sample(data, model, params, f, meta, fit_par)

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
            util.make_mcmc_rprs_txt(vals_mcmc, errs_lower_mcmc, errs_upper_mcmc, meta)

        if not meta.s30_fit_white and ('rp' in meta.labels) and meta.run_nested:
            plots.nested_rprs(vals_nested, errs_lower_nested, errs_upper_nested, meta)
            util.make_nested_rprs_txt(vals_nested, errs_lower_nested, errs_upper_nested, meta)

    print('Finished s30')

    return meta

def update_clips(clips_array):
    clips_old = clips_array[0]
    clips_new = clips_array[1]
    clips_new[0] + sum([i <= clips_new[0] for i in clips_old])
    clips_new_updated = [clips_new[ii] + sum([i <= clips_new[ii] for i in clips_old]) for ii in range(len(clips_new))]
    return np.concatenate((clips_old, clips_new_updated))


def labels_gen(params, meta, fit_par):
    nvisit = int(meta.nvisit)
    labels = []
    ii = 0
    for i in range(int(len(params) / nvisit)):
        if fit_par['fixed'][ii].lower() == "false":
            if str(fit_par['tied'][ii]) == "-1":
                labels.append(fit_par['parameter'][ii])
                ii = ii + 1
            else:
                for j in range(nvisit):
                    labels.append(fit_par['parameter'][ii] + str(j))
                    ii = ii + 1
        else:
            ii = ii + 1
    return labels
