import numpy as np
import sys, os, time
import unittest

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
    print(dirs)
    # saves times when these subdirectories were created. 
    # They always have the following form: 'run_YYYY-MM-DD_HH-MM-SS_eventlabel'
    dirs_bool = np.array(['/run_2' in i for i in dirs])
    dirs = dirs[dirs_bool]
    print(dirs)
    eventlabel_len = len(eventlabel)
    dirs_times = [i[-(eventlabel_len+20):-(eventlabel_len+1)] for i in dirs]
    print(dirs_times)
    # sort the times
    times_sorted = sn.sort_nicely(dirs_times)
    print(times_sorted)
    # most recent time
    recent_time = times_sorted[-1]
    print(recent_time)
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


class Testing(unittest.TestCase):
    def test_s00(self):
        print('test_path', test_path)

        pcf_path = test_path + '/run_files'

        #run s00
        meta = s00.run00(eventlabel, pcf_path)
        workdir = meta.workdir + '/'
        time.sleep(2)

        # run assertions for s00
        self.assertTrue(os.path.exists(workdir))
        self.assertTrue(os.path.exists(workdir+'/figs'))

    def test_s01(self):
        workdir, eventlabel = workdir_finder()
        time.sleep(2)

        #run s01
        meta = s01.run01(eventlabel, workdir)

        # run assertions for s01
        self.assertTrue(os.path.exists(workdir))
        self.assertTrue(os.path.exists(workdir+'/figs'))


if __name__ == '__main__':
    unittest.main()
