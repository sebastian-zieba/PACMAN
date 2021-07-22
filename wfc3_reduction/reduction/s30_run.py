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


def run30(eventlabel, workdir, meta=None):

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    #myfuncs = ['constant', 'upstream_downstream', 'ackbar', 'polynomial1', 'transit']
    myfuncs = ['constant', 'upstream_downstream', 'model_ramp', 'polynomial1', 'transit'] 
    #myfuncs = ['constant', 'upstream_downstream', 'model_ramp', 'polynomial2', 'transit']


    #defaults for command line flags
    verbose         = meta.run_verbose
    output          = meta.run_output
    show_plot       = meta.run_show_plot
    run_mcmc        = meta.run_mcmc
    nested          = meta.run_nested
    run_lsq         = meta.run_lsq
    plot_raw_data   = meta.run_plot_raw_data
    path            = workdir +"/extracted_lc/2021-07-22_06:58/"
    fit_white       = meta.run_fit_white
    divide_white    = meta.run_divide_white



    #reads in observation and fit parameters
    fit_par =   ascii.read(meta.workdir + "/fit_par.txt", Reader=ascii.CommentedHeader)

    files = ['/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-07-22_06-58-22_L-98-59_Hubble15856/extracted_lc/2021-07-22_06:58/lc_white.txt']#meta.run_files #glob.glob(os.path.join(path, "*"))
    if fit_white: files = glob.glob(white_file)

    meta.run_out_name = "fit_" + pythontime.strftime("%Y_%m_%d_%H:%M") + ".txt"

    for f in files:
        data = Data(f, meta, fit_par)
        model = Model(data, myfuncs)
        data, model, params = lsq_fit(fit_par, data, meta, model, myfuncs)

        """ind = model.resid/data.err > 10.
        print "num outliers", sum(ind)
        data.err[ind] = 1e12
        data, model = lsq_fit(fit_par, data, flags, model, myfuncs)"""
    
        ##rescale error bars so reduced chi-squared is one
        """data.err *= np.sqrt(model.chi2red)                                      
        data, model, params = lsq_fit(fit_par, data, flags, model, myfuncs)"""
        if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
        

        #FIXME : make this automatic!
        """outfile = open("white_systematics.txt", "w")
        for i in range(len(model.all_sys)): print(model.all_sys[i], file = outfile)
        outfile.close()"""
                            
        if meta.run_mcmc:
            output = mcmc_fit(data, model, params, f, meta, fit_par)

        if meta.run_nested:
            ##rescale error bars so reduced chi-squared is one
            data.err *= np.sqrt(model.chi2red)
            data, model, params = lsq_fit(fit_par, data, meta, model, myfuncs)
            if meta.run_verbose == True: print("rms, chi2red = ", model.rms, model.chi2red)
            output = nested_sample(data, model, params, f, meta, fit_par)

    return meta