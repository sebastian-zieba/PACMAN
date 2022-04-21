.. _stages:

Stages
============

PACMAN is separated into different stages, summarized here.  Click on the hyperlinks below to see the python script for each stage. For an example of each stage applied to real data, check out :ref:`example_introduction`.


- `Stage 00 <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s00_table.html>`_

  + read in  the _ima_ fits files in the data directory 
  + save important information from the header in a table
  + create a "work directory" to save all plots and files from the next stages


- `Stage 01: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s01_horizons.html>`_

  + download the location of HST during the observations to prepare for the barycentric correction


- `Stage 02: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s02_barycorr.html>`_

  + perform the barycentric correction to convert from MJD (which is in the header) to BJD


- `Stage 03: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s03_refspectra.html>`_

  + download or generate a stellar spectrum that matches the spectral type of the host star.  Options are: KURUCZ 1993, CASTELLI AND KURUCZ 2004, PHOENIX MODELS and a blackbody spectrum
  + multiply the stellar spectrum by the grism throughput to make a reference spectrum for the wavelength calibration 
  + saves the bandpass in pacman/ancil/bandpass


- `Stage 10: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s10_direct_images.html>`_

  + read in the direct images 
  + fit the position of the star in each direct image


- `Stage 20: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s20_extract.html>`_

  + extract the spectra using the optimal extraction routine from `Horne 1986 <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H>`_
  + save the white light curve into "lc_white.txt" and the spectroscopic light curve information into "lc_spec.txt".


- `Stage 21: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s21_bin_spectroscopic_lc.html>`_

  + bins the spectra in wavelength to create spectroscopic light curves (this step can be skipped if the user is only interested in the broadband light curve).


- `Stage 22: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s22_ld_inputmaker.html>`_

  + (optional, to be added) generate limb darkening parameters from a stellar model

- `Stage 30: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s30_run.html>`_

  + fit a model to the light curve(s)
