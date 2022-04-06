import numpy as np
import sys, os, time
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

    with capsys.disabled():
        print("PACMAN test:")

    pcf_path = test_path + '/run_files'

    #run s00
    meta = s00.run00(eventlabel, pcf_path)
    workdir = meta.workdir + '/'
    time.sleep(1.1)

    # run assertions for s00
    assert os.path.exists(workdir)
    assert os.path.exists(workdir+'/figs')

    filelist_path = workdir + '/filelist.txt'
    if os.path.exists(filelist_path):
        filelist = ascii.read(filelist_path)

    #print(t_mjd)
    ncols = len(filelist[0])
    nrows = len(filelist['t_mjd'])

    assert filelist['t_mjd'][0] == 56364.52973075
    assert (nrows, ncols) == (3, 9) 


@pytest.mark.dependency(depends=['test_s00'])
def test_s01(capsys):
    reload(s01)
    time.sleep(2)

    workdir, eventlabel = workdir_finder()

    #run s01
    meta = s01.run01(eventlabel, workdir)

    # run assertions for s01
    horizons_file = workdir+'/ancil/horizons/horizons_results_v0.txt' 

    assert os.path.exists(horizons_file)


@pytest.mark.dependency(depends=['test_s01'])
def test_s02(capsys):
    reload(s02)
    time.sleep(2)

    workdir, eventlabel = workdir_finder()

    #run s01
    meta = s02.run02(eventlabel, workdir)

    filelist_path = workdir + '/filelist.txt'
    if os.path.exists(filelist_path):
        filelist = ascii.read(filelist_path)

    # run assertions for s02
    print(filelist.colnames)
    assert ('t_bjd' in filelist.colnames)
    assert filelist['t_bjd'][0] == 2456365.030605767


@pytest.mark.dependency(depends=['test_s02'])
def test_s03(capsys):
    reload(s03)
    time.sleep(2)

    workdir, eventlabel = workdir_finder()

    #run s01
    meta = s03.run03(eventlabel, workdir)

    sm_file = workdir +  '/ancil/stellar_models/k93models/kp03_3500.fits'

    assert os.path.exists(sm_file)

    refspec_file = workdir + '/ancil/refspec/refspec.txt'

    assert os.path.exists(refspec_file)

    wvl_refspec, flux_refspec = np.loadtxt(refspec_file).T
    print(len(wvl_refspec))


    # run assertions for s02
    assert len(wvl_refspec) == 162
    assert flux_refspec[0] == 0

    #os.system("rm -r ./{0}".format(workdir))
    #return 0


