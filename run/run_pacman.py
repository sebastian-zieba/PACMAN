import sys, os, time
sys.path.append('/home/zieba/Desktop/Projects/Open_source/PACMAN/')

import pacman.reduction.s00_table as s00
import pacman.reduction.s01_horizons as s01
import pacman.reduction.s02_barycorr as s02
import pacman.reduction.s03_refspectra as s03
import pacman.reduction.s10_direct_images as s10
import pacman.reduction.s20_extract as s20
import pacman.reduction.s21_bin_spectroscopic_lc as s21
import pacman.reduction.s30_run as s30
from pacman.lib.update_meta import update_meta
from pacman.lib import sort_nicely as sn

run_s00 = False     # read in fits files and create filelist.txt
run_s01 = False     # download positions of HST during observations
run_s02 = False     # correct the MJD to BJD using the positions of HST
run_s03 = False      # download the stellar spectrum and create a reference spectrum with the bandpass of the grism
run_s10 = False     # determine the position of the source by looking at the direct image
run_s20 = False     # extract the spectra
run_s21 = False     # bin light curves
run_s30 = True     # fit models to the extracted light curve(s)

#eventlabel = 'KELT11_Hubble15926'
eventlabel = 'GJ1214_Hubble13021'
#workdir = '/home/zieba/Desktop/Projects/Open_source/PACMAN/run/run_2022-01-13_00-13-00_KELT11_Hubble15926/'

if run_s00:
    meta = s00.run00(eventlabel)
    runs = [i for i in os.listdir('.') if os.path.isdir(i)]
    runs = sn.sort_nicely(runs)
    workdir = runs[-1] + '/'

if not run_s00:
    runs = [i for i in os.listdir('.') if os.path.isdir(i)]
    runs = sn.sort_nicely(runs)
    workdir = runs[-1] + '/'

if run_s01:
    update_meta(eventlabel, workdir)
    meta = s01.run01(eventlabel, workdir)

if run_s02:
    update_meta(eventlabel, workdir)
    meta = s02.run02(eventlabel, workdir)

if run_s03:
    update_meta(eventlabel, workdir)
    meta = s03.run03(eventlabel, workdir)

if run_s10:
    update_meta(eventlabel, workdir)
    meta = s10.run10(eventlabel, workdir)

if run_s20:
    update_meta(eventlabel, workdir)
    meta = s20.run20(eventlabel, workdir)

if run_s21:
    update_meta(eventlabel, workdir)
    meta = s21.run21(eventlabel, workdir)

if run_s30:
    update_meta(eventlabel, workdir)
    meta = s30.run30(eventlabel, workdir)
