.. _before_running:

Before Running
================

Before we can run ``PACMAN`` we need to create a run directory.

1) **Decide on a "run directory"**

All runs will be saved in a run directory. For the analysis of the GJ1214 data, I have the following directory setup:

::

	user
	└── Desktop
		└── Projects
			└── Observations
				└── Hubble
					└── GJ1214_13021
						├── pacman_script.py
						├── fit_par.txt
						└── obs_par.pcf
		└── Data
			└── GJ1214_Hubble13021
                            ├── ibxy06d0q_ima.fits
                            ├── ...
                            └── ibxy07ryq_ima.fits

In this case `/user/Desktop/Projects/Observations/Hubble/GJ1214_13021` is my run directory.

2) **Copy the contents of run_files into your run directory**

The contents of run_files are three different files:

 - **pacman_script.py**

 - **fit_par.txt**

 - **obs_par.pcf**

These three files are `part of the package <https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files>`_ but can also be downloaded under this link: `Download here <https://downgit.github.io/#/home?url=https://github.com/sebastian-zieba/PACMAN/tree/master/src/pacman/data/run_files>`_.

3) **Set up the pcf**

The pacman control file, or :ref:`pcf page <pcf>`, stores all the parameters the user can tune for a particular run.
Before we start, paths to the data and the run directory (``data_dir`` and ``run_dir``) have to be set in the pcf.
Follow along to :ref:`Stage 00 <stage00>` to set these parameters.
