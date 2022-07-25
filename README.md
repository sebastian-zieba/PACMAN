![](https://github.com/sebastian-zieba/PACMAN/blob/master/docs/source/media/Pacman_V2.gif)

[![DOI](https://zenodo.org/badge/278655981.svg)](https://zenodo.org/badge/latestdoi/278655981)
[![Documentation Status](https://readthedocs.org/projects/pacmandocs/badge/?version=latest)](https://pacmandocs.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/sebastian-zieba/PACMAN/branch/master/graph/badge.svg?token=YGPOSJSH5Z)](https://codecov.io/gh/sebastian-zieba/PACMAN)


# PACMAN

**Welcome to PACMAN**
``PACMAN`` is a pipeline to reduce and analyze Hubble/Wide Field Camera 3 (WFC3) observations of transiting exoplanets. The pipeline runs end-to-end, beginning with a time series of 2D images and ending with a spectrum for the planet. ``PACMAN`` can easily fit multiple observations simultaneously.
                                                                                
The main steps in the pipeline are:                                             
                                                                                
- optimally extract spectra from the 'ima' data products provided by STScI      
- bin the spectra into user-specified wavelength bins and output the light curve(s)
- fit the light curves with a variety of astrophysical models (transit, eclipse, phase curve) and instrument systematic models (visit-long quadratic trends, orbit-long exponential trends)
- estimate uncertainties on the planet parameters with least-squares, MCMC, or nested sampling
                                                                                

## Installation and Documentation

Check out the [docs](https://pacmandocs.readthedocs.io/en/latest/) for the latest instructions for how to install and use ``PACMAN``.


## Having trouble?

Please take a look at the [issues page](https://github.com/sebastian-zieba/PACMAN/issues) to see if anyone else has encountered similar problems. If not, raise an issue and we will do our best to address it!
