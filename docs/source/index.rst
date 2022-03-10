.. PACMAN documentation master file, created by
   sphinx-quickstart on Mon Dec 20 11:05:56 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PACMAN's documentation!
==================================

.. image:: media/Pacman_V1.gif

**Welcome to the documentation for PACMAN**

``PACMAN`` is a pipeline to reduce and analyze HST data (G102 or G141).

It takes 'ima' data products created by STSci and outputs a white or spectroscopic light curves.
The user can then fit various functions (transit, eclipse, visit long quadratic trends,
orbital exponential trends, etc ...) to these light curves.

The pipeline is still under development! Use at your own risk.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   stages
   pcf
   directories
   models
   api
   contact

.. toctree::
   :maxdepth: 1
   :caption: Example GJ1214 b

   example_introduction
   download_data
   before_running
   stage00
   stage01
   stage02
   stage03
   stage10
   stage20
   stage21
   stage30



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
