"""This code computes the mean position of the direct image for each visit"""
from pathlib import Path

import numpy as np
from astropy.io import ascii, fits
from astropy.table import Table
from tqdm import tqdm

from .lib import util
from .lib import gaussfitter
from .lib import plots
from .lib import manageevent as me


def run10(eventlabel, workdir: Path, meta=None):
    """Opens the direct images to determine the position of the star on the detector.
    The positions are then saved in x and y physical pixel coordinates into a new txt file called xrefyref.txt.
    """
    print('Starting s10')
    if meta is None:
        meta = me.loadevent(workdir / f'WFC3_{eventlabel}_Meta_Save')

    #f = open(meta.workdir + '/xrefyref.txt', 'w')						#opens file to store positions of reference pixels
    table = Table() # creates table to store positions of reference pixels

    # load in more information into meta
    meta = util.ancil(meta, s10=True)

    t_bjd = meta.t_bjd_di
    iorbit_di = meta.iorbit_di
    ivisit_di = meta.ivisit_di

    pos1_all = np.zeros(len(meta.files_di))
    pos2_all = np.zeros(len(meta.files_di))

    files_sp = meta.files_sp     # spectra files
    sp = fits.open(files_sp[0])  # opens first spectral file

    #iterate over the direct images
    for i, file in enumerate(tqdm(meta.files_di, desc='Determining Source Positions for Direct Images', ascii=True)):
        ima = fits.open(file)

        #NEW: account for difference in CRPIX1 between direct image and spectroscopic images
        LTV1 = ima[1].header['LTV1'] + int(abs(ima[1].header['CRPIX1'] - sp[1].header['CRPIX1']))					#X offset to get into physical pixels + difference in CRPIX1 between direct and spectral image 
        LTV2 = ima[1].header['LTV2'] + int(abs(ima[1].header['CRPIX1'] - sp[1].header['CRPIX1']))					#Y offset to get to physical pixels + difference in CRPIX1 between direct and spectral image

        dat = ima[1].data[meta.di_rmin:meta.di_rmax, meta.di_cmin:meta.di_cmax]				#cuts out stamp around the target star
        err = ima[2].data[meta.di_rmin:meta.di_rmax, meta.di_cmin:meta.di_cmax]

        plots.image_quick(ima, i, meta)

        # If the guess for the cutout is outside of the dimensions of the dataset we will do the next iteration.
        # If this wouldn't be checked for, data AND err WILL BE EMPTY AND THE GAUSSFIT WILL TERMINATE THE STAGE.
        if meta.di_rmax > ima[1].data.shape[0] or meta.di_rmax < 0:
            print('\nYour guess for di_rmax is outside of the image.')
            continue
        if meta.di_cmax > ima[1].data.shape[1] or meta.di_cmax < 0:
            print('\nYour guess for di_cmax is outside of the image.')
            continue
        if meta.di_rmin > ima[1].data.shape[0] or meta.di_rmin < 0:
            print('\nYour guess for di_rmin is outside of the image.')
            continue
        if meta.di_cmin > ima[1].data.shape[1] or meta.di_cmin < 0:
            print('\nYour guess for di_cmin is outside of the image.')
            continue

        # run the fitter
        results = gaussfitter.gaussfit(dat, err)

        if meta.save_image_plot or meta.show_image_plot:
            plots.image(dat, ima, results, i, meta)

        # save positions into a file
        # TODO: convert this txt file to an astropy table
        # TODO: instead of listing the time and positions, the visit&orbit number and positions would be better
        pos1 = results[3]+meta.di_rmin-LTV1
        pos2 = results[2]+meta.di_cmin-LTV2
        pos1_all[i] = pos1
        pos2_all[i] = pos2
        #print(t_bjd[i], pos1, pos2, file=f)

        ima.close()
    #f.close()

    table['t_bjd']        = np.array(meta.t_bjd_di, dtype=np.float64)
    table['ivisit']       = np.array(meta.ivisit_di, dtype=np.int64)
    table['iorbit']       = np.array(meta.iorbit_di, dtype=np.int64)
    table['iorbit_cumul'] = np.array(meta.iorbit_di_cumulative, dtype=np.int64)
    table['pos1']	      = np.array(pos1_all, dtype=np.float64)
    table['pos2']	      = np.array(pos2_all, dtype=np.float64)
    ascii.write(table, meta.workdir / 'xrefyref.txt', format='rst', overwrite=True)

    util.di_reformat(meta)

    # Save results
    print('Saving Metadata')
    me.saveevent(meta, meta.workdir / f'WFC3_{meta.eventlabel}_Meta_Save', save=[])
    print('tmp: ', 2+4)
    print('Finished s10 \n')
    return meta
