![](https://github.com/sebastian-zieba/PACMAN/blob/master/docs/source/media/Pacman_V2.gif)

[![DOI](https://zenodo.org/badge/278655981.svg)](https://zenodo.org/badge/latestdoi/278655981)
[![status](https://joss.theoj.org/papers/73b451e3272d57cc63d1be78c097ec4d/status.svg)](https://joss.theoj.org/papers/73b451e3272d57cc63d1be78c097ec4d)
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


## Citing PACMAN
Please cite ``Zieba & Kreidberg (2022)`` if you use ``PACMAN`` in your research. 
The BibTeX entry for the paper is::

    @ARTICLE{2022JOSS....7.4838Z,
           author = {{Zieba}, Sebastian and {Kreidberg}, Laura},
            title = "{PACMAN: A pipeline to reduce and analyze Hubble Wide Field Camera 3 IR Grism data}",
          journal = {The Journal of Open Source Software},
         keywords = {astronomy, exoplanets, python, Python, spectroscopy, HST, photometry, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
             year = 2022,
            month = dec,
           volume = {7},
           number = {80},
              eid = {4838},
            pages = {4838},
              doi = {10.21105/joss.04838},
    archivePrefix = {arXiv},
           eprint = {2212.11421},
     primaryClass = {astro-ph.IM},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2022JOSS....7.4838Z},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }


## Issues and Contributing

If you are encountering issues with the code, please take a look at the [issues page](https://github.com/sebastian-zieba/PACMAN/issues) to see if anyone else has encountered similar problems (make sure you check both open and closed GitHub Issues!). 
If not, raise an issue and we will do our best to address it! 
When submitting an issues, please provide information like the current ``PACMAN`` version and the error message you are receiving in your terminal.
Also, have a look at our [FAQ page](https://pacmandocs.readthedocs.io/en/latest/faq.html) where we collect frequently asked questions. 
If you are thinking about contributing to ``PACMAN``, please also create an issue, so that we can make sure nobody else is already working on the same feature.
