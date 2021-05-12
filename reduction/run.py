import sys
sys.path.insert(0, './util')
sys.path.insert(0, './util/models')
import numpy as np
import glob
import os
from astropy.io import ascii
import getopt
import time as pythontime
from read_data import Data
from model import Model
from least_squares import lsq_fit
from mcmc import mcmc_fit
from nested import nested_sample
import yaml

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
    verbose         = False
    output          = False
    show_plot       = False
    run_mcmc        = False
    nested          = False
    run_lsq         = True
    plot_raw_data   = False
    path            = "/extracted_lc/2021-04-21_19:42/"
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
    obs_par_path = "config/obs_par.yaml"
    with open(obs_par_path, 'r') as file:
        obs_par = yaml.safe_load(file)
    fit_par =   ascii.read("config/fit_par.txt", Reader=ascii.CommentedHeader)

    files = ['/home/zieba/Desktop/Projects/Open_source/wfc3-pipeline/reduction/extracted_lc/2021-05-03_04:22/lc_white.txt']#glob.glob(os.path.join(path, "*"))
    if fit_white: files = glob.glob(white_file)

    flags['out-name'] = "fit_" + pythontime.strftime("%Y_%m_%d_%H:%M") + ".txt"

    for f in files:
        data = Data(f, obs_par, fit_par)
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
            output = mcmc_fit(data, model, params, f, obs_par, fit_par)     

        if flags['nested']:
            ##rescale error bars so reduced chi-squared is one
            data.err *= np.sqrt(model.chi2red)
            data, model, params = lsq_fit(fit_par, data, flags, model, myfuncs)
            if flags['verbose'] == True: print("rms, chi2red = ", model.rms, model.chi2red)
            output = nested_sample(data, model, params, f, obs_par, fit_par)


if __name__ == '__main__':
    main()
