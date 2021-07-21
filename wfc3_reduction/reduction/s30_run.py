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



def usage():
    cmd = sys.argv[0]
    sys.stderr.write('Usage: python %s OPTION\n\n' % os.path.basename(cmd))
    sys.stderr.write(
        'Allowed OPTION flags are:\n'
        '  --show-plot       displays fitted light curve plots\n'
        '  --run-mcmc        runs MCMC starting from least-squares parameters\n'
        '  --nested        runs dynesty\n'
        '  --plot-raw-data   plots raw light curve separated by visit\n'   
        '  --path PATH	     specifies PATH to light curves to be fit\n' 
        '  --fit-white FILE  fits the white light curve stored in FILE\n' 
        '  -v                prints fit diagnostic information\n'
        '  -o                saves fit output to file\n'
        '  -h                lists instructions for usage\n'
        '\n')
    sys.exit(1)


def main():


    eventlabel = 'L-98-59_Hubble15856'
    workdir = '/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-07-21_07-37-23_L-98-59_Hubble15856/'

    meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")

    #myfuncs = ['constant', 'upstream_downstream', 'ackbar', 'polynomial1', 'transit']
    myfuncs = ['constant', 'upstream_downstream', 'model_ramp', 'polynomial1', 'transit'] 
    #myfuncs = ['constant', 'upstream_downstream', 'model_ramp', 'polynomial2', 'transit'] 

    #significance above which to mask outliers
    #outlier_cut = 10.

    #parses command line input
    try: opts, args = \
            getopt.getopt(sys.argv[1:], 
                "hov", ["help", "show-plot", "run-mcmc", "nested", "plot-raw-data", 
                "plot-sys", "path=", "fit-white=", "divide-white"]
            )
    except getopt.GetoptError: usage()

    #defaults for command line flags
    verbose         = True
    output          = False
    show_plot       = True
    run_mcmc        = True
    nested          = False
    run_lsq         = True
    plot_raw_data   = True
    path            = "/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-07-21_07-37-23_L-98-59_Hubble15856/extracted_lc/2021-07-21_07:37"
    fit_white       = False
    divide_white    = False

    for o, a in opts:
        if o in ("-h", "--help"): usage()
        elif o == "-o": output = True
        elif o == "-v": verbose = True
        elif o == "--show-plot": show_plot = True
        elif o == "--run-mcmc": run_mcmc, run_lsq = True, False
        elif o == "--nested": nested, run_lsq = True, True
        elif o == "--run-lsq": run_lsq = True
        elif o == "--plot-raw-data": plot_raw_data = True
        elif o == "--path": path = a
        elif o == "--fit-white": fit_white, white_file = True, a
        elif o == "--divide-white": divide_white = True
        else: assert False, "unhandled option"

    flags = {'verbose': verbose, 'show-plot': show_plot, 
            'plot-raw-data': plot_raw_data, 'output': output, 
            'out-name': 'none.txt', 'run-lsq': run_lsq, 
            'run-mcmc': run_mcmc, 'nested': nested, 'divide-white': divide_white, 
            'fit-white': fit_white}		

    #reads in observation and fit parameters
    fit_par =   ascii.read(meta.workdir + "/fit_par.txt", Reader=ascii.CommentedHeader)

    files = ['/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/run/run_2021-07-21_07-37-23_L-98-59_Hubble15856/extracted_lc/2021-07-21_07:37/lc_white.txt']#glob.glob(os.path.join(path, "*"))
    if fit_white: files = glob.glob(white_file)

    flags['out-name'] = "fit_" + pythontime.strftime("%Y_%m_%d_%H:%M") + ".txt"

    for f in files:
        data = Data(f, meta, fit_par)
        model = Model(data, myfuncs)
        data, model, params = lsq_fit(fit_par, data, flags, model, myfuncs)

        """ind = model.resid/data.err > 10.
        print "num outliers", sum(ind)
        data.err[ind] = 1e12
        data, model = lsq_fit(fit_par, data, flags, model, myfuncs)"""
    
        ##rescale error bars so reduced chi-squared is one
        """data.err *= np.sqrt(model.chi2red)                                      
        data, model, params = lsq_fit(fit_par, data, flags, model, myfuncs)"""
        if flags['verbose'] == True: print("rms, chi2red = ", model.rms, model.chi2red)
        

        #FIXME : make this automatic!
        """outfile = open("white_systematics.txt", "w")
        for i in range(len(model.all_sys)): print(model.all_sys[i], file = outfile)
        outfile.close()"""
                            
        if flags['run-mcmc']:
            output = mcmc_fit(data, model, params, f, meta, fit_par)

        if flags['nested']:
            ##rescale error bars so reduced chi-squared is one
            data.err *= np.sqrt(model.chi2red)
            data, model, params = lsq_fit(fit_par, data, flags, model, myfuncs)
            if flags['verbose'] == True: print("rms, chi2red = ", model.rms, model.chi2red)
            output = nested_sample(data, model, params, f, meta, fit_par)


if __name__ == '__main__':
    main()
