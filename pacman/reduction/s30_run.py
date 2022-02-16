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
from ..lib import manageevent as me
from ..lib.read_data import Data
from ..lib.model import Model
from ..lib.least_squares import lsq_fit
from ..lib.mcmc import mcmc_fit
from ..lib.nested import nested_sample
from ..lib.formatter import ReturnParams
from ..lib import sort_nicely as sn
from ..lib import nice_fit_par
from ..lib import plots
from ..lib import util


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

    #read in the files (white or spectroscopic) which will be fitted
    if meta.s30_fit_white:
        print('White light curve fit will be performed')
        files = []
        if meta.s30_most_recent_s20:
            lst_dir = os.listdir(meta.workdir + "/extracted_lc/")
            lst_dir = sn.sort_nicely(lst_dir)
            white_dir = lst_dir[-1]
            files.append(meta.workdir + "/extracted_lc/" + white_dir + "/lc_white.txt")
        else:
            files.append(meta.s30_white_file)
    elif meta.s30_fit_spec:
        print('Spectroscopic light curve fit(s) will be performed')
        meta.wavelength_list = []
        if meta.s30_most_recent_s21:
            #TODO: This also includes the "bins12" when sorting. But i only want to sort by the time in the dir names. So just use the end of the dir names instead
            lst_dir = os.listdir(meta.workdir + "/extracted_sp/")
            lst_dir = sn.sort_nicely(lst_dir)
            spec_dir = lst_dir[-1]
            files = glob.glob(os.path.join(meta.workdir + "/extracted_sp/" + spec_dir, "*.txt"))  #
            files = sn.sort_nicely(files)
        else:
            files = glob.glob(os.path.join(meta.s30_spec_dir, "*.txt"))  #
            files = sn.sort_nicely(files)
    else:
        print('Neither s30_fit_white nor s30_fit_spec are True!')

    print('Identified file(s) for fitting:', files)

    if meta.run_verbose:
        vals = []
        errs = []
        idxs = []

    meta.chi2red_list = []

    for counter, f in enumerate(files):
        print('\n****** File: {0}/{1}'.format(counter+1, len(files)))
        meta.s30_file_counter = counter
        time.sleep(1.1) #sleep to prevent overwriting of data if saved in the same second
        meta.run_file = f
        meta.fittime = time.strftime('%Y-%m-%d_%H-%M-%S')

        if meta.run_clipiters == 0:
            print('\n')
            data = Data(f, meta, fit_par)
            model = Model(data, myfuncs)
            data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True) #not clipping
            meta.chi2red_list.append(model.chi2red)
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


        """ind = model.resid/data.err > 10.
        print "num outliers", sum(ind)
        data.err[ind] = 1e12
        data, model = lsq_fit(fit_par, data, flags, model, myfuncs)"""

        if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)

        #FIXME : make this automatic!
        """outfile = open("white_systematics.txt", "w")
        for i in range(len(model.all_sys)): print(model.all_sys[i], file = outfile)
        outfile.close()"""

        meta.labels = labels_gen(params, meta, fit_par)

        if meta.run_mcmc:
            time.sleep(1.1)
            if meta.rescale_uncert:
                ##rescale error bars so reduced chi-squared is one
                data.err *= np.sqrt(model.chi2red)
            data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
            if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
            mcmc_fit(data, model, params, f, meta, fit_par)

        if meta.run_nested:
            time.sleep(1.1)
            if meta.rescale_uncert:
                ##rescale error bars so reduced chi-squared is one
                data.err *= np.sqrt(model.chi2red)
            data, model, params, m = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
            if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
            nested_sample(data, model, params, f, meta, fit_par)

        if meta.run_verbose:
            val, err, idx = ReturnParams(m, data)
            vals.append(val)
            errs.append(err)
            idxs.append(idx)

    if meta.run_verbose:
        plots.params_vs_wvl(vals, errs, idxs, meta)

        if not meta.s30_fit_white:
            # Saves rprs and wvl as a txt file
            util.make_lsq_rprs_txt(vals, errs, idxs, meta)
            # Saves rprs vs wvl as a plot
            plots.lsq_rprs(rprs_vals_lsq, rprs_errs_lsq, meta)

        if not meta.s30_fit_white and meta.run_mcmc:
            plots.mcmc_rprs(meta)
            util.make_mcmc_rprs_txt(meta)

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
