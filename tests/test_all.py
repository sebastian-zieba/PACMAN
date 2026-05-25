import time
import os
import sys
from pathlib import Path
import shutil

import numpy as np
import pytest
from astropy.io import ascii, fits
from astroquery.mast import Observations
from astropy.table import Table
from importlib import reload
from photutils.datasets import make_noise_image

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
from pacman.lib import util


src_dir = Path.cwd() / 'src'
sys.path.insert(0, src_dir)

# HACK: Sets the encoding to utf-8 for the Github tests.
sys.stdout.reconfigure(encoding='utf-8')

test_path = Path(__file__).parent
pcf_path = test_path / "pacman_run_files"


def latest_stage(stage_num: str) -> Path:
    return util.find_latest_stage_run(
        test_path,
        f"stage{stage_num}",
        f"s{stage_num}_run_*",
    )


def delete_dir(dir_name: Path) -> None:
    if dir_name.exists():
        print(f"Old dir found and deleted: {dir_name}")
        shutil.rmtree(dir_name)


def delete_stage_dirs() -> None:
    for stage_dir in [
        "stage00", "stage01", "stage02", "stage03",
        "stage10", "stage20", "stage21", "stage30",
    ]:
        delete_dir(test_path / stage_dir)


def replace_pcf_value(path: Path, key: str, value: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    new_lines = []

    for line in lines:
        stripped = line.strip()
        if stripped.startswith(key):
            comment = ""
            if "#" in line:
                before, comment = line.split("#", 1)
                comment = " # " + comment.strip()
            new_lines.append(f"{key:<30} {value}{comment}")
        else:
            new_lines.append(line)

    path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")


def replace_fit_par_row(
    path: Path,
    parameter: str,
    fixed: str,
    tied: str,
    value: str,
    prior: str,
    p1: str,
    p2: str,
) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    new_lines = []

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("#") or not stripped:
            new_lines.append(line)
            continue

        parts = line.split()

        if parts[0] == parameter:
            new_lines.append(
                f"{parameter:<14}{fixed:<7}{tied:<6}"
                f"{value:<12}{prior:<7}{p1:<12}{p2:<12}".rstrip()
            )
        else:
            new_lines.append(line)

    path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")


@pytest.mark.run(order=1)
def test_sessionstart(capsys):
    """Called as the first test. It downloads the three HST files
    used in this test using astroquery.
    """
    test_dir = Path(__file__).parent
    delete_stage_dirs()

    # Delete old data dir
    data_dir = test_dir / 'data'
    mast_dir = test_dir / 'mastDownload'  # Specify root directory to be searched for .sav files.
    delete_dir(data_dir)
    delete_dir(mast_dir)

    # Create a data dir
    data_dir.mkdir(parents=True, exist_ok=True)

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
    for tree, _, fils in os.walk(mast_dir):
        filelist.extend([Path(tree) / fil for fil in fils if fil.endswith('.fits')])

    for fil in filelist:
        # HACK: On windows file will be not overwritten
        # so if it exists we remove it
        if (data_dir / fil.name).exists():
            os.remove(data_dir / fil.name)
        fil.rename(data_dir / fil.name)
    delete_dir(mast_dir)
    assert True


@pytest.mark.run(order=2)
def test_s00(capsys):
    """Reads in the downloaded HST files and creates the work directory
    and the filelist file.
    """
    reload(s00)
    pcf_path = test_path / 'pacman_run_files'

    # Run s00
    meta = s00.run00(pcf_path)
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

    # Run s01
    _ = s01.run01(pcf_path)
    s01_workdir = latest_stage("01")

    horizons_file = s01_workdir / 'ancil' / 'horizons' / 'horizons_results_v0.txt'

    # Run assertions
    assert horizons_file.exists()


def my_round(num):
    """Cutoff a decimal number after 2 decimal places without rounding.
    From: https://stackoverflow.com/questions/967661/python-truncate-after-a-hundreds
    """
    return float("%.2f" % (int(num*100)/float(100)))


@pytest.mark.run(order=4)
def test_horizons(capsys):
    """Check the shape of the HORIZONS file."""

    s01_workdir = latest_stage("01")
    horizons_file = s01_workdir / "ancil" / "horizons" / "horizons_results_v0.txt"

    start_data = '$$SOE'
    end_data = '$$EOE'

    # Read in whole table as an list of strings, one string per line
    with horizons_file.open('r') as ctable:
        wholetable = ctable.readlines()

    # Find start and end line
    i = 0
    # While end has not been found:
    while wholetable[i].find(end_data) == -1:
        # If start is found get the index of next line:
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

    # Run s02
    _ = s02.run02(pcf_path)
    s02_workdir = latest_stage("02")

    filelist_file = s02_workdir / 'filelist.txt'
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

    # Run s03
    _ = s03.run03(pcf_path)
    s03_workdir = latest_stage("03")

    sm_file = s03_workdir / 'ancil' / 'stellar_models' / 'k93models' / 'kp03_3500.fits'
    assert sm_file.exists()

    with fits.open(sm_file) as hdul:
        wvl = hdul[1].data["WAVELENGTH"] * 1e-10
        flux = hdul[1].data["g50"] * 1e-7 * 1e4 / 1e-10 / np.pi

    # Check if for the sm fits file the flux and wavelength is >= 0 everywhere
    assert np.all(wvl >= 0)
    assert np.all(flux >= 0)

    # Check the refspec_file
    refspec_file = s03_workdir / "ancil" / "refspec" / "refspec.txt"
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

    # Run s10
    _ = s10.run10(pcf_path)
    s10_workdir = latest_stage("10")


    xrefyref_file = s10_workdir / 'xrefyref.txt'
    assert xrefyref_file.exists()

    if xrefyref_file.exists():
        xrefyref = ascii.read(xrefyref_file)

    # Check if the direct image position was determined correctly
    assert np.isclose(xrefyref['pos1'][0], 513.57510, atol=1e-4)
    assert np.isclose(xrefyref['pos2'][0], 400.90239, atol=1e-4)


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

    yy, xx = np.mgrid[0:tshape[0], 0:tshape[1]]

    source_data = source['flux'][0] * np.exp(
        -0.5 * (
            ((xx - source['x_mean'][0]) / source['x_stddev'][0])**2
            + ((yy - source['y_mean'][0]) / source['y_stddev'][0])**2
        )
    )
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

    # Run s20
    meta = s20.run20(pcf_path)
    s20_workdir = latest_stage("20")

    s20_lc_spec_file = s20_workdir / "extracted_lc" / "lc_spec.txt"
    s20_lc_white_file = s20_workdir / "extracted_lc" / "lc_white.txt"
    s20_wvl_table_file = s20_workdir / "extracted_lc" / "wvl_table.dat"

    # Check if the files were created
    assert s20_lc_spec_file.exists()
    assert s20_lc_white_file.exists()
    assert s20_wvl_table_file.exists()

    s20_lc_spec = ascii.read(s20_lc_spec_file)
    s20_lc_white = ascii.read(s20_lc_white_file)

    # Check the amount of columns
    assert len(s20_lc_spec.colnames) == 10
    assert len(s20_lc_white.colnames) == 11

    s20_wvl_table = ascii.read(s20_wvl_table_file)
    assert s20_wvl_table.colnames == [
        "bin",
        "wavelength",
        "half_width",
        "lower_edge",
        "upper_edge",
    ]
    assert len(s20_wvl_table) == 1

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

    [f_opt_0, _, numoutliers] = optextr.optextr(
            spectrum, err, spec_box_0, var_box_0, Mnew,
            meta.nsmooth, meta.sig_cut, meta.save_optextr_plot, 0, 0, meta)

    assert np.round(np.sum(f_opt_0), 0) == np.round(np.sum(spectrum), 0) #optimal extraction flux should be the same as the total flux in the array
    assert numoutliers == 0  # We didnt introduce any outliers


@pytest.mark.run(order=29)
def test_s21(capsys):
    """Creates spectroscopic light curves."""
    reload(s21)
    time.sleep(1)

    # Run s21
    meta = s21.run21(pcf_path)
    s21_workdir = latest_stage("21")
    extracted_sp_dir_path = s21_workdir / "extracted_sp"

    s21_dir = s21_workdir / "extracted_sp"

    s21_wvl_table_file = s21_dir / "wvl_table.dat"
    assert s21_wvl_table_file.exists()

    s21_wvl_table = ascii.read(s21_wvl_table_file)

    assert s21_wvl_table.colnames == [
        "bin",
        "wavelength",
        "half_width",
        "lower_edge",
        "upper_edge",
    ]

    wvl_s21 = s21_wvl_table["wavelength"]

    # Check if the number of bins defined in the pcf is the same as
    # the number of wavelength bins saved into the wvl_table.dat file.
    assert meta.wvl_bins == len(wvl_s21)

    # Number of light curves should be the same as meta.wvl_bins
    extracted_sp_lcs_files = list(s21_dir.glob("speclc*.txt"))
    assert meta.wvl_bins == len(extracted_sp_lcs_files)

    # There should be 10 columns as for the /lc_spec.txt file which was generated after running s20.
    extracted_sp_lc_file_0 = sn.sort_nicely(extracted_sp_lcs_files)[0]
    extracted_sp_lc_0 = ascii.read(extracted_sp_lc_file_0)
    assert len(extracted_sp_lc_0.colnames) == 10


@pytest.mark.run(order=30)
def test_s30(capsys):
    """Fits spectroscopic and white light curves."""
    reload(s30)
    time.sleep(1)

    pcf_file = pcf_path / "obs_par.pcf"
    fit_par_file = pcf_path / "fit_par.txt"

    original_pcf = pcf_file.read_text(encoding="utf-8")
    original_fit_par = fit_par_file.read_text(encoding="utf-8")

    try:
        # ------------------------------------------------------------
        # 1) Spectroscopic fit
        # ------------------------------------------------------------
        replace_pcf_value(pcf_file, "s30_fit_white", "False")
        replace_pcf_value(pcf_file, "s30_fit_spec", "True")

        s21_workdir = latest_stage("21")
        extracted_sp_dir_path = s21_workdir / "extracted_sp"

        s21_dir = sorted([
            path for path in extracted_sp_dir_path.iterdir()
            if path.is_dir()
        ])[-1]

        spec_lc_files = sn.sort_nicely(list(s21_dir.glob("speclc*.txt")))

        s30.run30(pcf_path)

        s30_spec_workdir = util.find_latest_stage_run(
            test_path,
            "stage30/spec_lc",
            "s30_run_*",
        )

        assert s30_spec_workdir.exists()
        assert (s30_spec_workdir / "fit_spec").exists()
        assert (s30_spec_workdir / "fit_spec" / "fit_lc").exists()
        assert (s30_spec_workdir / "fit_spec" / "lsq_res").exists()

        # ------------------------------------------------------------
        # 2) White-light fit
        # ------------------------------------------------------------
        replace_pcf_value(pcf_file, "s30_fit_white", "True")
        replace_pcf_value(pcf_file, "s30_fit_spec", "False")

        replace_fit_par_row(
            fit_par_file,
            parameter="c",
            fixed="False",
            tied="-1",
            value="8.4",
            prior="U",
            p1="8.35",
            p2="8.45",
        )

        s30.run30(pcf_path)

        s30_white_workdir = util.find_latest_stage_run(
            test_path,
            "stage30/white_lc",
            "s30_run_*",
        )

        assert s30_white_workdir.exists()
        assert (s30_white_workdir / "fit_white").exists()
        assert (s30_white_workdir / "fit_white" / "fit_lc").exists()
        
    finally:
        pcf_file.write_text(original_pcf, encoding="utf-8")
        fit_par_file.write_text(original_fit_par, encoding="utf-8")

@pytest.mark.run(order=40)
def test_sessionfinish(capsys):
    """Called after whole test run finished. It will delete the created
    work directory and the downloaded HST files.
    """
    delete_stage_dirs()
    delete_dir(test_path / "data")
    delete_dir(test_path / "mastDownload")

    print('Deleted directories and files again.')
    assert True
