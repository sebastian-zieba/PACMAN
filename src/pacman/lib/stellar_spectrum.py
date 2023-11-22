import urllib.request
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from numpy.typing import ArrayLike

from .options import OPTIONS


def get_bb(user_teff):
    """Creates a blackbody spectrum for a given stellar effective
    temperature, Teff.

    Parameters
    ----------
    user_teff: float
        stellar effective temperature

    Returns
    -------
    wvl: numpy array
        wavelength np.linspace(0.1, 6, 1000) / 1e6
    flux: numpy array
        stellar flux in units of W/sr/m^3
    """
    wvls1 = np.linspace(0.1, 0.7, 100, endpoint=False) / 1e6
    wvls2 = np.linspace(0.7, 1.1, 400, endpoint=False) / 1e6
    wvls3 = np.linspace(1.1, 1.9, 400, endpoint=False) / 1e6
    wvls4 = np.linspace(1.9, 6.0, 100, endpoint=False) / 1e6
    wvls = np.concatenate((wvls1, wvls2, wvls3, wvls4))
    h = 6.62607004e-34  # m^2 kg/s
    c = 299792458.0  # m/s
    kb = 1.38064852e-23  # m^2 kg /s^2 K
    flux = (2.0*h*(c**2)/(wvls**5))*(1.0/( np.exp( (h*c)/(wvls*kb*user_teff) )- 1.0))   #SI units: W/sr/m^3
    return wvls, flux


def downloader(url: str) -> None:
    """This function downloads a file from the given url using
    urllib.request."""
    file_name = Path(url).name
    print(f'\t      + Downloading file {file_name} from {url}.')
    urllib.request.urlretrieve(url, file_name)


def find_nearest(array: ArrayLike, value: float):
    """Finds nearest element to a value in an array.

    Taken from https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def get_sm(meta, user_met, user_logg: float, user_teff: float):
    """Creates a Kurucz 1994, Castelli and Kurucz 2004 or Phoenix stellar
    spectrum for a given stellar effective temperature, metallicity and log g.

    Parameters
    ------------
    meta :
        a metadata instance.
    user_met : float
        stellar metallicity.
    user_logg : float
        stellar logg.
    user_teff : float
        stellar effective temperature.

    Returns
    ----------
    wvl : numpy array
        wavelength np.linspace(0.1, 6, 1000) / 1e6.
    flux : numpy array
        stellar flux in units of W/sr/m^3.

    Notes:
    ----------
    History:
        Written by Sebastian Zieba      December 2021
    """
    rooturl = 'https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/'
    sm = meta.sm
    if sm == 'k93models':
        label = 'k'
        possible_mets = np.array(
                [1.0, 0.5, 0.3, 0.2, 0.1, 0.0, -0.1,
                 -0.2, -0.3, -0.5, -1.0, -1.5, -2.0,
                 -2.5, -3.0, -3.5, -4.0, -4.5, -5.0])

    if sm == 'ck04models':
        label = 'ck'
        possible_mets = np.array([0.0, -0.5, -1.0, -1.5,
                                  -2.0, -2.5, +0.5, +0.2])

    if sm == 'phoenix':
        label = 'phoenix'
        possible_mets = np.array([0.0, -0.5, -1.0, -1.5, -2.0,
                                  -2.5, -3.0, -3.5, -4.0, +0.5, +0.3])

    chosen_met = find_nearest(possible_mets, user_met)
    if user_met not in possible_mets:
        print('Possible metallicities: {0}'.format(possible_mets))
        print('For input metallicity {}, closest metallicity is {}. \n'.format(user_met, chosen_met))
    elif user_met in possible_mets:
        print('Possible metallicities: {0}'.format(possible_mets))
        print('Using input metallicity of {}. \n'.format(user_met))

    # annoyingly, k93models and ck04models use the convention that 0 = +0 -> p00
    # but phoenix uses 0 = -0 -> m00
    if (sm == 'k93models') or (sm == 'ck04models'):
        if chosen_met < 0:
            l0 = 'm'
        else:
            l0 = 'p'
    elif sm == 'phoenix':
        if chosen_met <= 0:
            l0 = 'm'
        else:
            l0 = 'p'

    # e.g., M/H = 0.5 --> 05
    # e.g., M/H = 2.0 --> 20
    l1 = str(abs(chosen_met)).replace(".", "")

    met_url = f'{label}{l0}{l1}'  # e.g., ckp05
    full_url = f'{rooturl}{sm}/{met_url}'

    sm_dir_run = meta.workdir / 'ancil' / 'stellar_models'
    if not sm_dir_run.exists():
        sm_dir_run.mkdir(parents=True)

    sm_dir_run = sm_dir_run / sm
    if not sm_dir_run.exists():
        sm_dir_run.mkdir(parents=True)

    sm_dir_pkg = meta.pacmandir / 'data' / 'stellar_models' / sm
    if not sm_dir_pkg.exists():
        sm_dir_pkg.mkdir(parents=True)

    # check if a list of all downloadable files exists. If not create it.
    file_list_path = sm_dir_pkg / 'file_list.txt'
    if not file_list_path.exists():
        # inspired by code in:
        # https://stackoverflow.com/questions/40543200/want-to-get-all-links-in-a-webpage-using-urllib-request
        # https://github.com/nespinoza/limb-darkening/blob/master/get_lds.py

        html = str(urlopen(full_url).read())

        # hyperlinks in html have the following form, e.g.,:
        # <a href="kp00_5000.fits">kp00_5000.fits</a>
        # so we first look for <a
        # then take the label between "> and <\a>
        # finally check that the label includes .fits in the name
        all_files = []
        while True:
            idx1 = html.find('<a ')
            if (idx1 == -1):
                break
            else:
                html = html[idx1:]
                idx2 = html.find('\">')
                idx3 = html.find('</a>')
                filename = html[idx2 + 2:idx3]
                if '.fits' in filename:
                    all_files.append(html[idx2 + 2:idx3])
            html = html[idx3 + 4:]

        with file_list_path.open('w', encoding=OPTIONS["encoding"]) as f:
            for file in all_files:
                f.write(file + '\n')
    else:
        all_files = np.loadtxt(file_list_path, dtype=str).T

    # Determine the possible temperature from the filenames
    possible_teffs = np.array([float(i.split('_')[1].split('.')[0]) for i in all_files])
    chosen_teff = find_nearest(possible_teffs, user_teff)
    if user_teff not in possible_teffs:
        print(f'Possible effective temperatures: {possible_teffs}')
        print(f'For input effective temperature {user_teff}, closest temperature is {chosen_teff}.\n')
    elif user_teff in possible_teffs:
        print(f'Possible effective temperatures: {possible_teffs}')
        print(f'Using input effective temperature of {user_teff}.\n')

    #we now know the metallicity and temperature and can download the needed file (it wont download it again if it already exists)
    filename = f'{met_url}_{int(chosen_teff)}.fits'
    full_url = f'{full_url}/{filename}'

    # If the file wasnt downloaded yet, download it. Then move it into meta.workdir + 'ancil/stellar_models/{0}/'.format(sm)
    filepath = sm_dir_run / filename
    print(f'Was the stellar model fits file called {filename} already downloaded?:', filepath.exists(), '\n')
    if not filepath.exists():
        downloader(full_url)
        Path(filename).rename(filepath)

    # the files contains all available log g for a specific metallicity and temperature
    hdul = fits.open(filepath)
    cols = hdul[1].data.names

    # finding out which log g are available
    possible_loggs = []
    for i in cols:
        if 'g' in i:
            possible_loggs.append(float(i[1] + '.' + i[2]))

    possible_loggs = np.array(possible_loggs)
    chosen_logg = find_nearest(possible_loggs, user_logg)
    if user_logg not in possible_loggs:
        print(f'Possible logg: {possible_loggs}'.format(possible_loggs))
        print(f'For input logg {user_logg}, closest logg is {chosen_logg}.\n')
    elif user_logg in possible_loggs:
        print(f'Possible logg: {possible_loggs}')
        print(f'Using input logg of {user_logg}.\n')

    chosen_logg_name = f'g{str(chosen_logg)[0]}{str(chosen_logg)[2]}'
    wvl = hdul[1].data['WAVELENGTH']*1e-10
    flux = hdul[1].data[chosen_logg_name]*1e-7*1e4/1e-10/np.pi
    # flux = flux * wvl
    hdul.close()
    return wvl, flux

# def get_phoenix_spiderman(Teff, logg, MH):
#     PHOENIX_DIR = '/home/zieba/Documents/PHOENIX-ACES-AGSS-COND-2011_R10000FITS_Z-0.0'
#     ftemplate = 'lte{teff:05d}-{logg:4.2f}-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
#     filename = os.path.join(PHOENIX_DIR, ftemplate.format(teff=Teff, logg=logg, z=MH))
#
#     # changing to si, W / m^3 / str
#     flux, h = fits.getdata(filename, ext=0, header=True)
#
#     flux = flux * 1e-7 * 1e6 / (np.pi)
#
#     crval = h['CRVAL1']
#     cdelt = h['CDELT1']
#     ctype = h['CTYPE1']
#
#     if ctype == 'AWAV-LOG':
#         wvl = (np.exp(crval + cdelt * np.arange(0, len(flux)))) * 1e-10
#     else:
#         print('ctype is not log! It  is {}'.format(ctype))
#
#     return wvl, flux
#
#
# def gen_phoenix(wvls, Teff, MH, logg):
#     sm = 'phoenix'
#     # Generate Stellar spectrum
#     sp = S.Icat(sm, Teff, MH, logg)
#     #print(sp.waveunits.name)
#     #print(sp.fluxunits.name)
#     sp_wvl = sp.wave*1e-10 #m
#     sp_flux = sp.flux*1e-7*1e4/1e-10/np.pi
#     #x = sp_wvl
#     #y = sp_flux
#     #f = interp1d(x, y, kind='cubic')
#     return sp_wvl, sp_flux #wvls, f(wvls)
#
# def gen_kurutz(wvls, Teff, MH, logg):
#     sm = 'ck04models'
#     # Generate Stellar spectrum
#     sp = S.Icat(sm, Teff, MH, logg)
#     sp_wvl = sp.wave*1e-10 #m
#     sp_flux = sp.flux*1e-7*1e4/1e-10/np.pi
#     #x = sp_wvl
#     #y = sp_flux
#     #f = interp1d(x, y, kind='cubic')
#     return sp_wvl, sp_flux #wvls, f(wvls)
