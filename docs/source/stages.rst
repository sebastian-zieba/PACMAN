.. _stages:

Stages
============

PACMAN is separated into different stages.
Here is a quick summary of them.
You can find an example of every stage applied on real observations in :ref:`example_introduction`.
You can click on the Stage to inspect the main python script for it.



- `Stage 00: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s00_table.html>`_
  Every _ima_ fits file in the data directory is being read in and
  important information from the header are being saved into a table.
  This stage will create a "work directory" and all plots and files from the next stages
  will be saved there.


- `Stage 01: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s01_horizons.html>`_
  To prepare for the barycentric correction,
  we download the location of HST during the observations.


- `Stage 02: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s02_barycorr.html>`_
  We perform the barycentric correction to convert from MJD (which is in the header) to BJD.


- `Stage 03: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s03_refspectra.html>`_
  We make a reference spectrum for the wavelength calibration later on.
  For that, we first download or generate a stellar spectrum.
  Options are: KURUCZ 1993, CASTELLI AND KURUCZ 2004, PHOENIX MODELS and a blackbody spectrum.
  The product of the stellar spectrum and the bandpass of the Grism (G102 or G141)
  is the reference spectrum. The bandpass is saved in pacman/ancil/bandpass.


- `Stage 10: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s10_direct_images.html>`_
  The direct images are being read in
  and the position of the star is being determined.
  If there were several direct images taken at the beginning of an orbit,
  the user can decide if they want to use the median or latest direct image position.


- `Stage 20: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s20_extract.html>`_
  We extract the spectra.
  We use optimal extraction for that as described in `Horne 1986 <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H>`_.
  The white light curve information will be saved into "lc_white.txt" and the spectroscopic light curve information into "lc_spec.txt".
  If the user is only interested into fitting the white light curve, he or she can skip to Stage 30.


- `Stage 21: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s21_bin_spectroscopic_lc.html>`_
  We bin the spectra to create spectroscopic light curves (this step can be skipped if the user is only interested to fit the white light curve).


- `Stage 22: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s22_ld_inputmaker.html>`_
  This is a feature planned to be added and not implemented yet!
  Generate limb darkening parameters
  (this step can be skipped if the user does not want to calculate limb-darkening coefficients for a specified stellar model in order to fix them during fitting).


- `Stage 30: <https://pacmandocs.readthedocs.io/en/latest/_modules/pacman/reduction/s30_run.html>`_
  Fit a model to the white or spectroscopic light curve(s).
