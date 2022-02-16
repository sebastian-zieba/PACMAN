.. _stages:

Stages
============

PACMAN is separated into different stages. Here is a quick summary before we look into each of them into more detail by
applying PACMAN on real observations:

- Stage 00: Every observation in the data directory is being read in and important information from the header are put into a table and saved into a new file.

- Stage 01: To prepare for the barycentric correction, we download the location of HST during the observations.

- Stage 02: We perform the barycentric correction to convert from MJD to BJD.

- Stage 03: We load in a stellar spectrum. The product of the stellar spectrum and the bandpass of the Grism (G102 or G141) will be used as a reference spectrum in the wavelength calibration.

- Stage 10: The direct images are being read in and the position of the star is being determined.

- Stage 20: We extract the spectra.

- Stage 21: We bin the spectra to create spectroscopic light curves (this step can be skipped if the user is only interested to fit the white light curve).

- Stage 22: IN WORK! Generate limb darkening parameters (this step can be skipped if the user does not want to calculate limb-darkening coefficients for a specified stellar model in order to fix them during fitting).

- Stage 30: Fit a model to the white or spectroscopic light curve(s).
