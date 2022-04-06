import numpy as np
import sys, os, time, glob
import pytest
from astropy.io import ascii

sys.path.insert(0, '/home/zieba/Desktop/Projects/Open_source/PACMAN/')
#print(sys.path)

from pacman.lib import util

from pacman.reduction import s00_table as s00
from pacman.reduction import s01_horizons as s01
from pacman.reduction import s02_barycorr as s02
from pacman.reduction import s03_refspectra as s03
from pacman.reduction import s10_direct_images as s10
from pacman.reduction import s20_extract as s20
from pacman.reduction import s21_bin_spectroscopic_lc as s21
from pacman.reduction import s30_run as s30

from pacman.lib import sort_nicely as sn

from importlib import reload

test_path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'

eventlabel='GJ1214_13021'

def workdir_finder():
    eventlabel='GJ1214_13021'
    # list subdirectories in the run directory
    dirs = np.array([f.path for f in os.scandir(test_path) if f.is_dir()])
    #print(dirs)
    # saves times when these subdirectories were created. 
    # They always have the following form: 'run_YYYY-MM-DD_HH-MM-SS_eventlabel'
    dirs_bool = np.array(['/run_2' in i for i in dirs])
    dirs = dirs[dirs_bool]
    #print(dirs)
    eventlabel_len = len(eventlabel)
    dirs_times = [i[-(eventlabel_len+20):-(eventlabel_len+1)] for i in dirs]
    #print(dirs_times)
    # sort the times
    times_sorted = sn.sort_nicely(dirs_times)
    #print(times_sorted)
    # most recent time
    recent_time = times_sorted[-1]
    #print(recent_time)
    # find the directory with that most recent time
    idx = 0
    for i in range(len(dirs)):
        if dirs[i][-(eventlabel_len+20):-(eventlabel_len+1)] == recent_time:
            idx = i
    workdir = dirs[idx]
    #save the eventlabel which is in the directory name too
    print('workdir: ', workdir)
    print('eventlabel: ', eventlabel)
    return (workdir, eventlabel)


@pytest.mark.dependency()
def test_s00(capsys):
    reload(s00)
    #print('test_path', test_path)
    #print(os.system("pwd"))

    pcf_path = test_path + '/run_files'

    #run s00
    meta = s00.run00(eventlabel, pcf_path)
    workdir = meta.workdir + '/'
    time.sleep(1)

    # run assertions
    assert os.path.exists(workdir)
    assert os.path.exists(workdir+'/figs')

    filelist_file = workdir + '/filelist.txt'
    assert os.path.exists(filelist_file)
    filelist = ascii.read(filelist_file)

    ncols = len(filelist[0])
    nrows = len(filelist['t_mjd'])

    assert filelist['t_mjd'][0] == 56364.52973075
    assert (nrows, ncols) == (3, 9) 


@pytest.mark.dependency(depends=['test_s00'])
def test_s01(capsys):
    reload(s01)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    #run s01
    meta = s01.run01(eventlabel, workdir)

    horizons_file = workdir+'/ancil/horizons/horizons_results_v0.txt' 

    # run assertions
    assert os.path.exists(horizons_file)


@pytest.mark.dependency(depends=['test_s01'])
def test_s02(capsys):
    reload(s02)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    #run s02
    meta = s02.run02(eventlabel, workdir)

    filelist_file = workdir + '/filelist.txt'
    assert os.path.exists(filelist_file)
    filelist = ascii.read(filelist_file)

    # Check if the barycentric correction was correctly performed
    assert ('t_bjd' in filelist.colnames)
    assert filelist['t_bjd'][0] == 2456365.030605767


@pytest.mark.dependency(depends=['test_s02'])
def test_s03(capsys):
    reload(s03)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    #run s03
    meta = s03.run03(eventlabel, workdir)

    sm_file = workdir +  '/ancil/stellar_models/k93models/kp03_3500.fits'

    assert os.path.exists(sm_file)

    refspec_file = workdir + '/ancil/refspec/refspec.txt'

    assert os.path.exists(refspec_file)

    wvl_refspec, flux_refspec = np.loadtxt(refspec_file).T

    # Check if the refspec was correctly created
    assert len(wvl_refspec) == 162
    assert flux_refspec[0] == 0


@pytest.mark.dependency(depends=['test_s03'])
def test_s10(capsys):
    reload(s10)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    #run s10
    meta = s10.run10(eventlabel, workdir)

    xrefyref_file = workdir + '/xrefyref.txt'

    assert os.path.exists(xrefyref_file)

    if os.path.exists(xrefyref_file):
        xrefyref = ascii.read(xrefyref_file)

    # Check if the direct image position was determined correctly
    assert np.round(xrefyref['pos1'][0], 5) == 513.57510
    assert np.round(xrefyref['pos2'][0], 5) == 400.90239


@pytest.mark.dependency(depends=['test_s10'])
def test_s20(capsys):
    reload(s20)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    #run s20
    meta = s20.run20(eventlabel, workdir)

    extracted_lc_dir_path = workdir + '/extracted_lc'

    s20_dir = np.array([f.path for f in os.scandir(extracted_lc_dir_path) if f.is_dir()])[0]

    s20_lc_spec_file = s20_dir + '/lc_spec.txt'
    s20_lc_white_file = s20_dir + '/lc_white.txt'

    #Check if the files were created
    assert os.path.exists(s20_lc_spec_file)
    assert os.path.exists(s20_lc_white_file)

    s20_lc_spec = ascii.read(s20_lc_spec_file)
    s20_lc_white = ascii.read(s20_lc_white_file)

    #Check the amount of columns
    assert len(s20_lc_spec.colnames) == 10
    assert len(s20_lc_white.colnames) == 11


@pytest.mark.dependency(depends=['test_s20'])
def test_s21(capsys):
    reload(s21)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    #run s21
    meta = s21.run21(eventlabel, workdir)

    extracted_sp_dir_path = workdir + '/extracted_sp'

    s21_dir = np.array([f.path for f in os.scandir(extracted_sp_dir_path) if f.is_dir()])[0]

    s21_wvl_table_file = s21_dir + '/wvl_table.dat'
    assert os.path.exists(s21_wvl_table_file)
    s21_wvl_table = ascii.read(s21_wvl_table_file)

    wvl_s21 = s21_wvl_table['wavelengths']

    #Check if the number of bins defined in the pcf is the same as 
    #the number of wavelength bins saved into the wvl_table.dat file.
    assert meta.wvl_bins == len(wvl_s21)

    #Number of light curves should be the same as meta.wvl_bins
    extracted_sp_lcs_files = glob.glob(os.path.join(s21_dir, "*.txt"))
    assert meta.wvl_bins == len(extracted_sp_lcs_files)

    #There should be 10 columns as for the /lc_spec.txt file which was generated after running s20.
    extracted_sp_lc_file_0 = sn.sort_nicely(extracted_sp_lcs_files)[0]
    extracted_sp_lc_0 = ascii.read(extracted_sp_lc_file_0)
    assert len(extracted_sp_lc_0.colnames) == 10


@pytest.mark.dependency(depends=['test_s21'])
def test_s30(capsys):
    reload(s30)
    time.sleep(1)

    workdir, eventlabel = workdir_finder()

    #run s30
    meta = s30.run30(eventlabel, workdir)

    workdir_dirs = np.array([f.path for f in os.scandir(workdir) if f.is_dir()])  
    fit_dirs = workdir_dirs[np.array(['fit_' in i for i in workdir_dirs])]
    fit_dir = fit_dirs[0]
    assert os.path.exists(fit_dir)
    #os.system("rm -r ./{0}".format(workdir))
    #return 0



