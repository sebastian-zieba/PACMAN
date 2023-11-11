import time
import os
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import ascii, fits
from astroquery.mast import Observations
from astropy.table import Table
from importlib import reload
from photutils.datasets import make_noise_image, make_gaussian_sources_image

from pacman import s00_table as s00
from pacman import s01_horizons as s01
from pacman import s02_barycorr as s02
from pacman import s03_refspectra as s03
from pacman import s10_direct_images as s10
from pacman import s20_extract as s20
from pacman import s21_bin_spectroscopic_lc as s21
from pacman import s30_run as s30
from pacman.lib import sort_nicely as sn
from pacman.lib.suntimecorr import getcoords as getcoords
from pacman.lib.gaussfitter import gaussfit as gaussfit
from pacman.lib import optextr

src_dir = Path.cwd() / 'src'
sys.path.insert(0, src_dir)

test_path = Path(__file__).parent
eventlabel = 'GJ1214_13021'


def workdir_finder():
    """Finds the latest work directory created.
    After running the Stage 00 test,
    we want to base all following tests (for s01, s02, ...) on the
    workdirectory created when running s00.
    """
    eventlabel = 'GJ1214_13021'
    # List subdirectories in the run directory
    dirs = np.array([path for path in test_path.iterdir() if path.is_dir()])

    # Saves times when these subdirectories were created.
    # They always have the following form: 'run_YYYY-MM-DD_HH-MM-SS_eventlabel'
    dirs_bool = np.array(['run_2' in dir.name for dir in dirs])
    dirs = dirs[dirs_bool]

    eventlabel_len = len(eventlabel)
    dirs_times = [i.name[-(eventlabel_len+20):-(eventlabel_len+1)] for i in dirs]

    # Sort the times
    times_sorted = sn.sort_nicely(dirs_times)

    # Most recent time
    recent_time = times_sorted[-1]

    # Find the directory with that most recent time
    idx = 0
    for i in range(len(dirs)):
        if dirs[i].name[-(eventlabel_len+20):-(eventlabel_len+1)] == recent_time:
            idx = i
    workdir = dirs[idx]
    # Save the eventlabel which is in the directory name too
    print('workdir: ', workdir)
    print('eventlabel: ', eventlabel)
    return workdir, eventlabel


def delete_dir(dir_name):
    if os.path.exists(dir_name):
        print('Old dir found and deleted')
        os.system("rm -r {0}".format(dir_name))


@pytest.mark.run(order=1)
def test_sessionstart(capsys):
    """Called as the first test. It downloads the three HST files
    used in this test using astroquery.
    """
    test_dir = Path(__file__).parent

    eventlabel = 'GJ1214_13021'
    dirs = np.array([path for path in test_path.iterdir() if path.is_dir()])
    dirs_bool = np.array(['run_2' in dir.name for dir in dirs])
    dirs = dirs[dirs_bool]
    for diri in dirs:
        delete_dir(diri)

    # Delete old data dir
    data_dir = test_dir / 'data'
    mast_dir = test_dir / 'mastDownload'  # Specify root directory to be searched for .sav files.
    delete_dir(data_dir)
    delete_dir(mast_dir)

    # Create a data dir
    data_dir.mkdir(parents=True)

    # Search for the HST data
    proposal_obs = Observations.query_criteria(proposal_id=13021,
                                               instrument_name='WFC3/IR', project='HST')
    data_products = Observations.get_product_list(proposal_obs)

    # NOTE: Just download these six files
    select = ['ibxy07p9q', 'ibxy07paq', 'ibxy07pbq',
              'ibxy07pcq', 'ibxy07pdq', 'ibxy07pfq']

    data_products_select = []
    for j in select:
        data_products_select.append((data_products['obs_id'] == j).data)
    data_products_new = data_products[np.any(data_products_select, axis=0)]
    data_products_ima = data_products_new[data_products_new['productSubGroupDescription'] == 'IMA']

    # NOTE: Download the three files
    Observations.download_products(data_products_ima, mrp_only=False,
                                   download_dir=str(test_dir))

    filelist = []
    for tree, fol, fils in os.walk(mast_dir):
        filelist.extend([Path(tree) / fil for fil in fils if fil.endswith('.fits')])

    for fil in filelist:
        fil.rename(data_dir / fil.name)
    os.system(f"rm -r {mast_dir}")
    assert True


@pytest.mark.run(order=2)
def test_s00(capsys):
    """Reads in the downloaded HST files and creates the work directory
    and the filelist file.
    """
    reload(s00)
    pcf_path = test_path / 'run_files'

    # Run s00
    meta = s00.run00(eventlabel, pcf_path)
    time.sleep(1)

    # run assertions
    assert meta.workdir.exists()
    assert (meta.workdir / 'figs').exists()

    filelist_file = meta.workdir / 'filelist.txt'
    assert filelist_file.exists()
    filelist = ascii.read(filelist_file)

    ncols = len(filelist[0])
    nrows = len(filelist['t_mjd'])

    assert np.round(filelist['t_mjd'][0], 4) == 56364.5297
    assert (nrows, ncols) == (6, 10)


@pytest.mark.run(order=3)
def test_s01(capsys):
    """Downloads the HORIZONS file."""
    reload(s01)
    time.sleep(1)
    workdir, eventlabel = workdir_finder()

    # Run s01
    meta = s01.run01(eventlabel, workdir)
    horizons_file = workdir / 'ancil' / 'horizons' / 'horizons_results_v0.txt'

    # run assertions
    assert horizons_file.exists()


def my_round(num):
    """Cutoff a decimal number after 2 decimal places without rounding.
    From: https://stackoverflow.com/questions/967661/python-truncate-after-a-hundreds
    """
    return float("%.2f" % (int(num*100)/float(100)))


@pytest.mark.run(order=4)
def test_horizons(capsys):
    """Check the shape of the HORIZONS file."""
    workdir, eventlabel = workdir_finder()

    horizons_file = os.path.join(workdir,'ancil','horizons','horizons_results_v0.txt' )

    start_data = '$$SOE'
    end_data = '$$EOE'

    # Read in whole table as an list of strings, one string per line
    ctable = open(horizons_file, 'r')
    wholetable = ctable.readlines()
    ctable.close()

    # Find start and end line
    i = 0
    # while end has not been found:
    while wholetable[i].find(end_data) == -1:
        # if start is found get the index of next line:
        if wholetable[i].find(start_data) != -1:
            start = i + 1
        i += 1

    # Chop table
    data = wholetable[start:i - 2]

    # Extract values:
    x, y, z, time = getcoords(data)

    # Checking shape
    assert len(x) == 27

    # Checking first and last values
    assert np.all(np.array([my_round(x[0]), my_round(y[0]), my_round(z[0])]) == np.array([-147684997.27, 16573698.09, 7180590.09]))
    assert np.all(np.array([my_round(x[-1]), my_round(y[-1]), my_round(z[-1])])== np.array([-147721652.49, 16371575.31, 7082911.34]))


@pytest.mark.run(order=10)
def test_s02(capsys):
    """Performs the barycentric correction."""
    reload(s02)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    # Run s02
    meta = s02.run02(eventlabel, workdir)

    filelist_file = workdir / 'filelist.txt'
    assert filelist_file.exists()
    filelist = ascii.read(filelist_file)

    # Check if the barycentric correction was correctly performed
    assert ('t_bjd' in filelist.colnames)
    assert np.round(filelist['t_bjd'][0], 4) == 2456365.0306


@pytest.mark.run(order=15)
def test_s03(capsys):
    """Downloads the stellar spectrum and multiplies it with the bandpass."""
    reload(s03)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    # Run s03
    meta = s03.run03(eventlabel, workdir)
    sm_file = workdir / 'ancil' / 'stellar_models' / 'k93models' / 'kp03_3500.fits'
    assert sm_file.exists()

    hdul = fits.open(sm_file)
    wvl = hdul[1].data['WAVELENGTH']*1e-10
    flux = hdul[1].data['g50']*1e-7*1e4/1e-10/np.pi

    # Check if for the sm fits file the flux and wavelength is >= 0 everywhere
    assert np.all(wvl >= 0)
    assert np.all(flux >= 0)

    # Check the refspec_file
    refspec_file = workdir / 'ancil' / 'refspec' / 'refspec.txt'
    assert refspec_file.exists()

    wvl_refspec, flux_refspec = np.loadtxt(refspec_file).T

    # Check if the refspec was correctly created
    assert len(wvl_refspec) == 7751
    # Check if for the refspec file the flux and wavelength is >= 0 everywhere
    assert np.all(wvl_refspec >= 0)
    assert np.all(flux_refspec >= 0)


@pytest.mark.run(order=16)
def test_s10(capsys):
    """Determines the position of the direct image."""
    reload(s10)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    # Run s10
    meta = s10.run10(eventlabel, workdir)

    xrefyref_file = workdir / 'xrefyref.txt'
    assert xrefyref_file.exists()

    if xrefyref_file.exists():
        xrefyref = ascii.read(xrefyref_file)

    # Check if the direct image position was determined correctly
    assert np.round(xrefyref['pos1'][0], 5) == 513.57510
    assert np.round(xrefyref['pos2'][0], 5) == 400.90239


@pytest.mark.run(order=26)
def test_sim_source(capsys):
    """Determines the position of the simulated direct image."""
    sigma_psf = 1.0
    source = Table()

    source_xpos = np.random.randint(10, 54)
    source_ypos = np.random.randint(10, 54)

    source['flux'] = [5000]
    source['x_mean'] = [source_xpos]
    source['y_mean'] = [source_ypos]
    source['x_stddev'] = sigma_psf * np.ones(1)
    source['y_stddev'] = source['x_stddev']
    source['theta'] = [0]
    source['id'] = [1]
    tshape = (64, 64)

    source_data = make_gaussian_sources_image(tshape, source)
    noise_data = make_noise_image(tshape, distribution='gaussian', mean=80.,
                                  stddev=5., seed=123)

    image = (source_data + noise_data)
    results = gaussfit(image, noise_data)

    assert source_xpos-1 <= results[2] <= source_xpos+1
    assert source_ypos-1 <= results[3] <= source_ypos+1


@pytest.mark.run(order=27)
def test_s20(capsys):
    """The extraction step. Extracts flux as a function
    of wavelength and time.
    """
    reload(s20)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    # Run s20
    meta = s20.run20(eventlabel, workdir)

    extracted_lc_dir_path = workdir / 'extracted_lc'

    s20_dir = np.array([path for path in extracted_lc_dir_path.iterdir() if path.is_dir()])[0]
    s20_lc_spec_file = s20_dir / 'lc_spec.txt'
    s20_lc_white_file = s20_dir / 'lc_white.txt'

    # Check if the files were created
    assert s20_lc_spec_file.exists()
    assert s20_lc_white_file.exists()

    s20_lc_spec = ascii.read(s20_lc_spec_file)
    s20_lc_white = ascii.read(s20_lc_white_file)

    # Check the amount of columns
    assert len(s20_lc_spec.colnames) == 10
    assert len(s20_lc_white.colnames) == 11

    # test_optextr
    spectrum = np.ones((20, 9))

    for i in range(len(spectrum)):
        for j in range(len(spectrum[0])):
            if 4 < i < 8:
                if 1 < j < 7:
                    spectrum[i, j] = 10

    err = np.ones((20, 9))*0.01
    Mnew = np.ones((20, 9))
    spec_box_0 = 15 * 10
    var_box_0 = 1

    [f_opt_0, var_opt_0, numoutliers] = optextr.optextr(spectrum, err, spec_box_0, var_box_0, Mnew, meta.nsmooth, meta.sig_cut, meta.save_optextr_plot, 0, 0, meta)

    assert np.round(np.sum(f_opt_0), 0) == np.round(np.sum(spectrum), 0) #optimal extraction flux should be the same as the total flux in the array 
    assert numoutliers == 0  # We didnt introduce any outliers


@pytest.mark.run(order=29)
def test_s21(capsys):
    """Creates spectroscopic light curves."""
    reload(s21)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    # Run s21
    meta = s21.run21(eventlabel, workdir)

    extracted_sp_dir_path = workdir / 'extracted_sp'

    s21_dir = np.array([path for path in extracted_sp_dir_path.iterdir()
                        if path.is_dir()])[0]
    s21_wvl_table_file = s21_dir / 'wvl_table.dat'
    assert s21_wvl_table_file.exists()
    s21_wvl_table = ascii.read(s21_wvl_table_file)

    wvl_s21 = s21_wvl_table['wavelengths']

    # Check if the number of bins defined in the pcf is the same as
    # the number of wavelength bins saved into the wvl_table.dat file.
    assert meta.wvl_bins == len(wvl_s21)

    # Number of light curves should be the same as meta.wvl_bins
    extracted_sp_lcs_files = list(s21_dir.glob("*.txt"))
    assert meta.wvl_bins == len(extracted_sp_lcs_files)

    # There should be 10 columns as for the /lc_spec.txt file which was generated after running s20.
    extracted_sp_lc_file_0 = sn.sort_nicely(extracted_sp_lcs_files)[0]
    extracted_sp_lc_0 = ascii.read(extracted_sp_lc_file_0)
    assert len(extracted_sp_lc_0.colnames) == 10


@pytest.mark.run(order=30)
def test_s30(capsys):
    """Fits spectroscopic light curves."""
    reload(s30)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    # Run s30
    meta = s30.run30(eventlabel, workdir)

    dirs = np.array([path for path in workdir.iterdir() if path.is_dir()])
    dirs_bool = np.array(['fit_' in dir.name for dir in dirs])
    fit_dirs = dirs[dirs_bool]
    fit_dir = fit_dirs[0]
    assert fit_dir.exists()

    meta.s30_fit_white = True
    meta.s30_most_recent_s20 = True

    s30.run30(eventlabel, workdir, meta=meta)

    dirs = np.array([path for path in workdir.iterdir() if path.is_dir()])
    dirs_bool = np.array(['fit_' in dir.name for dir in dirs])

    print('dirs_bool: ', dirs_bool)
    assert True


@pytest.mark.run(order=40)
def test_sessionfinish(capsys):
    """Called after whole test run finished. It will delete the created
    work directory and the downloaded HST files.
    """
    workdir, _ = workdir_finder()
    test_dir = Path(__file__).parent
    data_dir = test_dir / 'data'
    os.system(f"rm -r {data_dir}")
    os.system(f"rm -r {workdir}")

    print('Deleted directories and files again.')
    assert True
