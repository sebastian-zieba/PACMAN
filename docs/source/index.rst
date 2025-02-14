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


Citing PACMAN
'''''''''''''

If you use ``PACMAN`` in your research, please cite `Zieba & Kreidberg (2022) <https://ui.adsabs.harvard.edu/abs/2022JOSS....7.4838Z/abstract>`_. The BibTeX entry for the paper is:


.. code-block::

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
   before_running
   astroquery_visits
   stage00
   stage01
   stage02
   stage03
   stage10
   stage20
   stage21
   stage30
   outputs

.. toctree::
   :maxdepth: 1
   :caption: Others

   download_data


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
