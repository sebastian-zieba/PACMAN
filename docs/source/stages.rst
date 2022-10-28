.. _stages:

Stages
============

PACMAN is separated into different stages, summarized here.  Click on the hyperlinks below to see the python script for each stage. For an example of each stage applied to real data, check out :ref:`example_introduction`.


- Stage 00

  + read in  the _ima_ fits files in the data directory 
  + save important information from the header in a table
  + create a "work directory" to save all plots and files from the next stages


- Stage 01

  + download the location of HST during the observations to prepare for the barycentric correction


- Stage 02

  + perform the barycentric correction to convert from MJD (which is in the header) to BJD


- Stage 03

  + download or generate a stellar spectrum that matches the spectral type of the host star.  Options are: KURUCZ 1993, CASTELLI AND KURUCZ 2004, PHOENIX MODELS and a blackbody spectrum
  + multiply the stellar spectrum by the grism throughput (which is saved in `PACMAN/src/data/bandpass <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/bandpass>`_) to make a reference spectrum for the wavelength calibration
  + saves the reference spectrum into the work directory


- Stage 10

  + read in the direct images 
  + fit the position of the star in each direct image


- Stage 20

  + extract the spectra using the optimal extraction routine from `Horne 1986 <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H>`_
  + save the white light curve into "lc_white.txt" and the spectroscopic light curve information into "lc_spec.txt".


- Stage 21

  + bins the spectra in wavelength to create spectroscopic light curves (this step can be skipped if the user is only interested in the broadband light curve).


- Stage 30

  + fit a model to the light curve(s). You can find the models currently available `here <https://pacmandocs.readthedocs.io/en/latest/models.html>`_.
  + the user can also download and use precalculated limb darkening parameters in this stage. For more information see the `fix-ld <https://pacmandocs.readthedocs.io/en/latest/pcf.html#fix-ld>`_ parameter in the pcf.
