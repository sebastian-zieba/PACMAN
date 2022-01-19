.. _quickstart:

Quickstart
============

1. Installation and requirements
------------------------------------

See installation.rst


2. Download data
------------------------------------

ima files!


3. Setup the PACMAN control file (.pcf)
------------------------------------

See pcf.rst

4. Run PACMAN
------------------------------------

PACMAN is separated into different stages. Here a quick summary before we look into each of them ińto more detail by
applying PACMAN on real observations:

- Stage 00: Every observation in the data directory is being read in an import information from the header are saved
into a new file.
- Stage 01: To prepare for the barycentric correction, we download the location of HST during the observations
- Stage 02: We perform the barycentric correction
- Stage 03: We load in a stellar spectrum. The product of the stellar spectrum and the bandpass of the Grism
(G102 or G141) will be used as a reference spectrum in the wavelength calibration.
- Stage 10: The direct images are being read in and the position of the star is being determined.
- Stage 20: We extract the spectra
- Stage 21: TODO! Generate limb darkening parameters
- Stage 22: Bin the spectra to create light curves
- Stage 30: Fit a model to the light curves

Stage 00
:::::::::::::::::::::::::::::::::::::::::


