import numpy as np
from exotic_ld import StellarLimbDarkening
import os
from astropy.table import QTable
from tqdm import tqdm
from astropy.io import ascii


def get_ld(meta):
    # Stellar models: 1D or 3D grid.
    ld_model = '1D'

    # Path to the installed data.
    ld_data_path = '/home/zieba/Downloads/exotic-ld_data'

    sld = StellarLimbDarkening(meta.MH, meta.Teff, meta.logg, ld_model, ld_data_path)

    # Instrument mode.
    if meta.grism == 'G102':
        mode = 'HST_WFC3_G102'
    elif meta.grism == 'G141':
        mode = 'HST_WFC3_G141'

    meta.wvl_bins = int(meta.wvl_bins)
    wave_edges = np.linspace(meta.wvl_min, meta.wvl_max, meta.wvl_bins+1)*1e4

    all_lds = []

    ld_dir_run = meta.workdir + 'ancil/limb_darkening/'
    if not os.path.exists(ld_dir_run):
        os.mkdir(ld_dir_run)


    if meta.ld_model == 2:
        table_ld = QTable(names=('wave_lower', 'wave_mid', 'wave_upper', 'u1', 'u2'))
    elif meta.ld_model == 4:
        table_ld = QTable(names=('wave_lower', 'wave_mid', 'wave_upper', 'c1', 'c2', 'c3', 'c4'))

    for i in tqdm(range(len(wave_edges)-1), desc='-- Calculating limb-darkening parameteres', leave=True, position=0):
        # Start and end of wavelength interval [angstroms].
        wave_mid = (wave_edges[i] + wave_edges[i+1])/2.
        wavelength_range = [wave_edges[i], wave_edges[i+1]]

        if meta.ld_model == 2:
            u1, u2 = sld.compute_quadratic_ld_coeffs(wavelength_range, mode)
            all_lds.append((u1, u2))
            table_ld.add_row([wave_edges[i], wave_mid, wave_edges[i+1], u1, u2])
        elif meta.ld_model == 4:
            c1, c2, c3, c4 = sld.compute_4_parameter_non_linear_ld_coeffs(wavelength_range, mode)
            all_lds.append((c1, c2, c3, c4))
            table_ld.add_row([wave_edges[i], wave_mid, wave_edges[i+1], c1, c2, c3, c4])

    ld_path = ld_dir_run + 'ld_file.txt'

    ascii.write(table_ld, ld_path, format='rst', overwrite=True)

    return ld_path