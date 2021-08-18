import sys
sys.path.append('../..')
sys.path.append('/home/zieba/Desktop/Projects/Open_source/wfc-pipeline/')
sys.path.insert(0, '../lib')
import numpy as np
import glob
import os
from astropy.io import ascii
import getopt
import time as pythontime
from ..lib import manageevent as me
from ..lib.read_data import Data
from ..lib.model import Model
from ..lib.least_squares import lsq_fit
from ..lib.mcmc import mcmc_fit
from ..lib.nested import nested_sample
import time
import shutil


def run30(eventlabel, workdir, meta=None):

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")


    # Create directories for Stage 3 processing
    datetime = time.strftime('%Y-%m-%d_%H-%M-%S')
    meta.fitdir = '/fit_' + datetime + '_' + meta.eventlabel
    if not os.path.exists(meta.workdir + meta.fitdir):
        os.makedirs(meta.workdir + meta.fitdir)

    # Copy ecf
    shutil.copy(meta.workdir + "/obs_par.ecf", meta.workdir + meta.fitdir)
    shutil.copy(meta.workdir + "/fit_par.txt", meta.workdir + meta.fitdir)


    myfuncs = meta.run_myfuncs

    #reads in observation and fit parameters
    fit_par =   ascii.read(meta.workdir + "/fit_par.txt", Reader=ascii.CommentedHeader)
    if meta.run_fit_white:
        files = meta.run_files #
    else:
        files = glob.glob(os.path.join(meta.run_files[0], "*.txt"))  #
    print(files)

    meta.run_out_name = "fit_" + pythontime.strftime("%Y_%m_%d_%H:%M") + ".txt"

    for f in files:

        meta.fittime = time.strftime('%Y-%m-%d_%H-%M-%S')

        if meta.run_clipiters == 0:
            print('\n')
            data = Data(f, meta, fit_par)
            model = Model(data, myfuncs)
            data, model, params = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
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
                    data, model, params = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
                else:
                    data, model, params, clip_idx = lsq_fit(fit_par, data, meta, model, myfuncs)
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
                            
        if meta.run_mcmc:
            ##rescale error bars so reduced chi-squared is one
            data.err *= np.sqrt(model.chi2red)
            data, model, params = lsq_fit(fit_par, data, meta, model, myfuncs, noclip=True)
            if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
            output = mcmc_fit(data, model, params, f, meta, fit_par)

        if meta.run_nested:
            ##rescale error bars so reduced chi-squared is one
            data.err *= np.sqrt(model.chi2red)
            data, model, params = lsq_fit(fit_par, data, meta, model, myfuncs)
            if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
            output = nested_sample(data, model, params, f, meta, fit_par)

    return meta

def update_clips(clips_array):
    clips_old = clips_array[0]
    clips_new = clips_array[1]

    clips_new[0] + sum([i <= clips_new[0] for i in clips_old])
    clips_new_updated = [clips_new[ii] + sum([i <= clips_new[ii] for i in clips_old]) for ii in range(len(clips_new))]

    return np.concatenate((clips_old, clips_new_updated))