.. PACMAN documentation master file, created by
   sphinx-quickstart on Mon Dec 20 11:05:56 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PACMAN's documentation!
==================================

.. image:: media/Pacman_V2.gif

**Welcome to the documentation for PACMAN**

``PACMAN`` is a pipeline to reduce and analyze Hubble/Wide Field Camera 3 (WFC3) observations of transiting exoplanets. The pipeline runs end-to-end, beginning with a time series of 2D images and ending with a spectrum for the planet. ``PACMAN`` can easily fit multiple observations simultaneously.

The main steps in the pipeline are:

- optimally extract spectra from the 'ima' data products provided by STScI 
- bin the spectra into user-specified wavelength bins and output the light curve(s)
- fit the light curves with a variety of astrophysical models (transit, eclipse, phase curve) and instrument systematic models (visit-long quadratic trends, orbit-long exponential trends) 
- estimate uncertainties on the planet parameters with least-squares, MCMC, or nested sampling

For a more detailed roadmap, see `Stages <https://pacmandocs.readthedocs.io/en/latest/stages.html>`_.

``PACMAN`` is being developed on GitHub. If you find a bug or want a new feature, please `raise an issue <https://github.com/sebastian-zieba/PACMAN/issues>`_!


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   stages
   directories
   pcf
   models
   api
   faq
   contact

.. toctree::
   :maxdepth: 1
   :caption: Quickstart

   example_introduction
   astroquery_visits
   before_running
   stage00
   stage01
   stage02
   stage03
   stage10
   stage20
   stage21
   stage30

.. toctree::
   :maxdepth: 1
   :caption: Others

   download_data


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
