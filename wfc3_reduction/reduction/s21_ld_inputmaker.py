import numpy as np
import yaml
from ..lib import manageevent as me
import os


def run21(eventlabel, workdir, meta=None):

    if meta == None:
        meta = me.loadevent(workdir + '/WFC3_' + eventlabel + "_Meta_Save")


    if meta.grism == 'G141':
        grism = 'g141_throughput.txt'
    elif meta.grism == 'G102':
        grism = 'g102_throughput.txt'

    Teff, logg, MH = meta.Teff, meta.logg, meta.MH


    for wvl_bins in meta.wvl_bins:
        #what bins do you want?
        #wave_bins = np.linspace(1.125, 1.65, 22)*1e4
        wave_bins = np.linspace(meta.wvl_min, meta.wvl_max, wvl_bins)*1e4
        print(wave_bins)

        dirname = meta.workdir + "/extracted_sp/"
        if not os.path.exists(dirname): os.makedirs(dirname)


        f = open(meta.workdir + "/extracted_sp/" + 'ld_inputfile_bins{}.txt'.format(wvl_bins), 'w')

        for i in range(wvl_bins-1):
            params = [i, Teff, logg, MH, '-1', grism, 'P100', wave_bins[i], wave_bins[i+1]]
            params_bin = '\t'.join(map(str,params))
            #print(params_bin)
            print(params_bin, file = f)
        f.close()

    print('Finished s21 \n')

    return meta