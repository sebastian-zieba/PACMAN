from pathlib import Path

import pacman.s00_table as s00
import pacman.s01_horizons as s01
import pacman.s02_barycorr as s02
import pacman.s03_refspectra as s03
import pacman.s10_direct_images as s10
import pacman.s20_extract as s20
import pacman.s21_bin_spectroscopic_lc as s21
import pacman.s30_run as s30

pcf_path = Path(__file__).parent

if __name__ == "__main__":
    meta = s00.run00(pcf_path=pcf_path) # reads in fits files and creates filelist.txt

    # meta = s01.run01(pcf_path=pcf_path) # downloads positions of HST during observations

    # meta = s02.run02(pcf_path=pcf_path) # corrects the MJD to BJD using the positions of HST

    # meta = s03.run03(pcf_path=pcf_path) # downloads the stellar spectrum and creates a reference spectrum with the bandpass of the grism

    # meta = s10.run10(pcf_path=pcf_path) # determines the position of the source by looking at the direct image

    # meta = s20.run20(pcf_path=pcf_path) # extracts the spectra

    # meta = s21.run21(pcf_path=pcf_path) # bins light curves

    # meta = s30.run30(pcf_path=pcf_path) # fits models to the extracted light curve(s)
