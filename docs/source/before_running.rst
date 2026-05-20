.. _before_running:

Before Running
==============

.. topic:: Summary

    - Decide on a data directory
    - Decide on a run directory
    - Create a ``pacman_run_files`` directory inside the run directory
    - Copy ``run_pacman.py``, ``fit_par.txt``, and ``obs_par.pcf`` into ``pacman_run_files``
    - Enter the paths to the run directory and data directory into the pcf


General workflow
-----------------

PACMAN automatically detects the most recent run directory from the previous stage.

For example:
- Stage 01 reads the latest ``stage00/s00_run_*``
- Stage 02 reads the latest ``stage01/s01_run_*``
- Stage 03 reads the latest ``stage02/s02_run_*``

and so on.

Each stage creates a new timestamped run directory while preserving previous runs.

Decide on a data directory
--------------------------

Create a directory where you would like to store the FITS files that will be reduced and analyzed.

In the next part of this Quickstart, we will
`download these FITS files <https://pacmandocs.readthedocs.io/en/latest/astroquery_visits.html>`_
and move them into this directory.


Decide on a run directory
-------------------------

Before we can run ``PACMAN``, we need to create a run directory.

All PACMAN outputs will be stored inside this directory. For the analysis of the
GJ 1214 data, we will use the following directory structure:

::

    user
    └── Desktop
        └── Projects
            └── Observations
                └── Hubble
                    └── GJ1214_13021  (="run directory")
                        └── pacman_run_files
                            ├── run_pacman.py
                            ├── fit_par.txt
                            └── obs_par.pcf
        └── Data
            └── GJ1214_Hubble13021  (="data directory")
                ├── ibxy06d0q_ima.fits
                ├── ...
                └── ibxy07ryq_ima.fits

In this example,

::

    /user/Desktop/Projects/Observations/Hubble/GJ1214_13021

is the run directory.


Copy the ``pacman_run_files`` templates
---------------------------------------

The PACMAN source code contains a directory called ``pacman_run_files``. It contains the following template files:

- ``run_pacman.py``
- ``fit_par.txt``
- ``obs_par.pcf``

These files are included in the PACMAN package and can also be downloaded here:

- `GitHub directory <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/pacman_run_files>`_
- `Direct download <https://downgit.github.io/#/home?url=https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/pacman_run_files>`_

Copy these files into the ``pacman_run_files`` directory inside the run directory.

You should then have the following structure:

::

    GJ1214_13021
    └── pacman_run_files
        ├── run_pacman.py
        ├── fit_par.txt
        └── obs_par.pcf


Set up the pcf
--------------

The PACMAN control file, or :ref:`pcf page <pcf>`, stores all parameters that the user can tune for a particular run.

Before starting the analysis, the paths to the data directory and run directory
(``datadir`` and ``rundir``) must be set in the pcf (to be clear, in the ``obs_par.pcf`` file that is in e.g., ``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021``).

Navigate into ``pacman_run_files`` (in e.g., ``/home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021``), open ``obs_par.pcf``, and modify the entries
``datadir`` and ``rundir`` with the absolute paths to the data directory and run directory, respectively.

For example:

| ``rundir   /home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021``
| ``datadir  /home/zieba/Desktop/Data/GJ1214_Hubble13021``


At the start of every stage, the current ``obs_par.pcf`` and ``fit_par.txt``
from ``pacman_run_files`` are copied into the new stage run directory.

The copied ``obs_par.pcf`` is then used to update the metadata for that stage.

This ensures every stage preserves the exact settings used for that run.



Run PACMAN
----------

(Don't do this yet. We will run PACMAN after downloading the data in the next part of this Quickstart.)
Navigate into the ``pacman_run_files`` directory and run:

.. code-block:: console

    python run_pacman.py

PACMAN will automatically create the stage directories (``stage00``,
``stage01``, etc.) inside the run directory as needed.


`In the next step of the quickstart <https://pacmandocs.readthedocs.io/en/latest/astroquery_visits.html>`_,
we will download the FITS files needed for the analysis using the Python package ``astroquery``.