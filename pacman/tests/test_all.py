import numpy as np
import sys, os, time
import pytest

sys.path.insert(0, '/home/zieba/Desktop/Projects/Open_source/PACMAN/')

from importlib import reload

from ..lib import util

from ..reduction import s00_table as s00
from ..reduction import s01_horizons as s01
from ..reduction import s02_barycorr as s02
from ..reduction import s03_refspectra as s03
from ..reduction import s10_direct_images as s10
from ..reduction import s20_extract as s20
from ..reduction import s21_bin_spectroscopic_lc as s21
from ..reduction import s30_run as s30

from ..lib import sort_nicely as sn

#from ..lib.util import pathdirectory

#@pytest.fixture(scope="session")


def test_s00(capsys):
    #print(os.system("pwd"))

    with capsys.disabled():
        print("PACMAN test:")

    # explicitly define meta variables to be able to run pathdirectory fn locally
    eventlabel='GJ1214_13021'
    #meta.topdir='../tests'
    pcf_path='./run_files/'

    #run s00
    meta = s00.run00(eventlabel, pcf_path)
    workdir = meta.workdir + '/'
    time.sleep(1.1)

    # run assertions for s00
    assert os.path.exists(workdir)
    assert os.path.exists(workdir+'/figs')


def test_s01(capsys):
    #print(os.system("pwd"))

    with capsys.disabled():
        print("PACMAN test:")

    # explicitly define meta variables to be able to run pathdirectory fn locally
    eventlabel='GJ1214_13021'
    #meta.topdir='../tests'
    pcf_path='./run_files/'


    # list subdirectories in the run directory
    dirs = np.array([f.path for f in os.scandir('.') if f.is_dir()])
    # saves times when these subdirectories were created. 
    # They always have the following form: 'run_YYYY-MM-DD_HH-MM-SS_eventlabel'
    dirs_bool = np.array([i[:7] == './run_2' for i in dirs])
    dirs = dirs[dirs_bool]
    print(dirs)
    dirs_times = [i[6:25] for i in dirs]
    # sort the times
    times_sorted = sn.sort_nicely(dirs_times)
    # most recent time
    recent_time = times_sorted[-1]
    # find the directory with that most recent time
    idx = 0
    for i in range(len(dirs)):
        if dirs[i][6:25] == recent_time:
            idx = i
    workdir = dirs[idx]
    #save the eventlabel which is in the directory name too
    eventlabel = workdir[26:]
    workdir = workdir + '/'
    print('workdir: ', workdir)
    print('eventlabel: ', eventlabel)


    #run s01
    meta = s01.run01(eventlabel, workdir)

    # run assertions for s01
    assert os.path.exists(workdir)
    assert os.path.exists(workdir+'/figs')


    #
    #meta = s02.run02(eventlabel, workdir)
    #meta = s03.run03(eventlabel, workdir)
    #meta = s10.run10(eventlabel, workdir)
    #meta = s20.run20(eventlabel, workdir)
    #meta = s21.run21(eventlabel, workdir)


    # remove temporary files
    #os.system("rm -r ./{0}".format(workdir))
    #os.system("rm -r data/JWST-Sim/NIRCam/Stage4/*")
    #os.system("rm -r data/JWST-Sim/NIRCam/Stage5/*")
