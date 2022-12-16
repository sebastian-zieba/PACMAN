.. _before_running:

Before Running
================

.. topic:: Summary
    - Decide on a data directory
    - Decide on a run directory
    - Copy pacman_script.py, fit_par.txt and obs_par.pcf to this run directory
    - Enter path to run directory and data directory into the pcf


Decide on a "data directory"
-----------------------------------
Create a new directory where you would like to store the fits files which we will reduce and analyze.
In the next part of this Quickstart, you will `download these fits files <https://pacmandocs.readthedocs.io/en/latest/astroquery_visits.html>`_ and move them into this directory.


Decide on a "run directory"
-----------------------------------

Before we can run ``PACMAN`` we need to create a run directory.

All runs will be saved in a run directory. For the analysis of the GJ1214 data, I have the following directory setup:

::

	user
	└── Desktop
		└── Projects
			└── Observations
				└── Hubble
					└── GJ1214_13021 (="run directory")
						├── pacman_script.py
						├── fit_par.txt
						└── obs_par.pcf
		└── Data
			└── GJ1214_Hubble13021 (="data directory")
                            ├── ibxy06d0q_ima.fits
                            ├── ...
                            └── ibxy07ryq_ima.fits

In this case, `/user/Desktop/Projects/Observations/Hubble/GJ1214_13021` is my run directory.


Copy run_files contents into your run directory
------------------------------------------------------------

The source code of PACMAN contains a directory called "run_files". In there are three different files:

 - **pacman_script.py**

 - **fit_par.txt**

 - **obs_par.pcf**

These three files are `part of the package <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files>`_ but can also be downloaded under this link: `Download here <https://downgit.github.io/#/home?url=https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files>`_.

This files have to be copied over to the run directory. You should then have the following:

::

	GJ1214_13021 (="run directory")
	├── pacman_script.py
	├── fit_par.txt
	└── obs_par.pcf


Set up the pcf
---------------------------------------

The pacman control file, or :ref:`pcf page <pcf>`, stores all the parameters the user can tune for a particular run.
Before we start, paths to the data and the run directory (``data_dir`` and ``run_dir``) have to be set in the pcf.

Navigate to the run directory, open the pcf and modify the entries ``data_dir`` and ``run_dir`` with the absolute paths to the data directory and run directory, respectively.
E.g.,

| ``rundir   /home/zieba/Desktop/Projects/Observations/Hubble/GJ1214_13021``
| ``datadir   /home/zieba/Desktop/Data/GJ1214_Hubble13021``


`In the next step of the quickstart <https://pacmandocs.readthedocs.io/en/latest/astroquery_visits.html>`_,
we will download the data (fits files) needed for our analysis using the python package ``astroquery``.


